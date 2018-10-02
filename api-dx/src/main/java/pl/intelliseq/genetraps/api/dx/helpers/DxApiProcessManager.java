package pl.intelliseq.genetraps.api.dx.helpers;

import com.dnanexus.*;
import com.dnanexus.exceptions.TagsException;
import com.dnanexus.exceptions.WrongNumberOfFilesException;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.core.env.Environment;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class DxApiProcessManager {

    @Autowired
    Environment env;

    public DXJob runUrlFetch(String inputUrl, String sampleid, String... tags) {
        var input = new HashMap<>();

        boolean leftTag = Arrays.asList(tags).contains("left");
        boolean rightTag = Arrays.asList(tags).contains("right");
        if(leftTag == true && rightTag == true) {
            throw new TagsException("both 'left' and 'right' tags cannot be used at same time");
        }
        if(leftTag == false && rightTag == false) {
            throw new TagsException("neither 'left' nor 'right' tag was used");
        }

        input.put("url", inputUrl);
        input.put("tags", tags);

        return DXApp.getInstance(env.getProperty("dx-url-fetch-ap-id"))
                .newRun()
                .setProject(DXProject.getInstance(env.getProperty("dx-project")))
                .setFolder(String.format("/samples/%s/rawdata", sampleid))
                .setInput(input)
                .run();
    }

    public void runMkDir(Integer sampleid) {
        runMkDir(sampleid.toString());
    }

    public void runMkDir(String sampleid) {
        DXContainer.getInstance(env.getProperty("dx-project")).newFolder("/samples/" + sampleid, true);
    }

    public DXJob runFastqc(String fileId) {
        var input = new HashMap<>();
        DXFile file = DXFile.getInstance(fileId);
        input.put("fastq_file", file);


        return getAppletFromName("iseq_fastqc")
                .newRun()
                .setProject(DXProject.getInstance(env.getProperty("dx-project")))
                .setInput(input)
                .setFolder(file.describe().getFolder().replace("rawdata", "fastqc/" + file.describe().getName()))
                .run();
    }

//    test me not in use
//    public DXJob runBwa(String fastq_fileId_1, String fastq_fileId_2) {
//        var input = new HashMap<>();
//        DXFile fastqFile1 = DXFile.getInstance(fastq_fileId_1);
//        DXFile fastqFile2 = DXFile.getInstance(fastq_fileId_2);
//        DXFile ref = DXFile.getInstance(env.getProperty("dx-reference"));
//
//        input.put("fastq_file_1", fastqFile1);
//        input.put("fastq_file_2", fastqFile2);
//        input.put("reference", ref);
//
//
//        return getAppletFromName("iseq_bwa")
//                .newRun()
//                .setProject(DXProject.getInstance(env.getProperty("dx-project")))
//                .setInput(input)
//                .setFolder(fastqFile1.describe().getFolder().replace("rawdata", "bwa"))
//                .run();
//    }

    public DXJob runBwa(Integer sampleid) {
        var input = new HashMap<>();
        List fileListLeft = DXSearch.findDataObjects().inFolder(DXContainer.getInstance(env.getProperty("dx-project")), "/samples/" + sampleid.toString() + "/rawdata").withTag("left").execute().asList();
        if(fileListLeft.size() != 1) {
            throw new WrongNumberOfFilesException("expected: 1 with 'left' tag ; found: " + fileListLeft.size());
        }
        DXFile fastqFile1 = (DXFile) fileListLeft.get(0);

        List fileListRight = DXSearch.findDataObjects().inFolder(DXContainer.getInstance(env.getProperty("dx-project")), "/samples/" + sampleid.toString() + "/rawdata").withTag("right").execute().asList();
        if(fileListRight.size() != 1) {
            throw new WrongNumberOfFilesException("expected: 1 with 'right' tag ; found: " + fileListRight.size());
        }
        DXFile fastqFile2 = (DXFile) fileListRight.get(0);
        DXFile ref = DXFile.getInstance(env.getProperty("dx-reference"));

        input.put("fastq_file_1", fastqFile1);
        input.put("fastq_file_2", fastqFile2);
        input.put("reference", ref);


        return getAppletFromName("iseq_bwa")
                .newRun()
                .setProject(DXProject.getInstance(env.getProperty("dx-project")))
                .setInput(input)
                .setFolder(fastqFile1.describe().getFolder().replace("rawdata", "bwa"))
                .run();
    }

    public DXJob runGatkHC(Integer sampleid, String interval) {
        var input = new HashMap<>();
        DXFile bam = (DXFile) DXSearch.findDataObjects().nameMatchesRegexp(".*\\.bam$")
                .inFolder(DXContainer.getInstance(env.getProperty("dx-project")), "/samples/" + sampleid.toString() + "/bwa")
                .execute().asList().get(0);
        DXFile ref = DXFile.getInstance(env.getProperty("dx-reference"));

        String[] vcfsIds = env.getProperty("dx-vcfsgz", String[].class);
        DXFile[] vcfs = new DXFile[vcfsIds.length];
        for(int i = 0; i < vcfs.length; i++) {
            vcfs[i] = DXFile.getInstance(vcfsIds[i]);
        }

        input.put("vcfs", vcfs);
        input.put("bam", bam);
        input.put("reference", ref);
        if(interval != null)
            input.put("interval", interval);


        return getAppletFromName("iseq_gatk_haplotype_caller")
                .newRun()
                .setProject(DXProject.getInstance(env.getProperty("dx-project")))
                .setInput(input)
                .setFolder(bam.describe().getFolder().replace("bwa", "gatkhc"))
                .run();
    }

    public DXApplet getAppletFromName(String appletName) {
        return (DXApplet) DXSearch.findDataObjects().nameMatchesExactly(appletName).execute().asList().get(0);
    }

    public JsonNode JSONDescribe(String id) {
        return DXJSON.safeTreeToValue(
                new DXHTTPRequest(DXEnvironment.create()).request("/" + id + "/" + "describe",
                        new ObjectMapper().createObjectNode(), DXHTTPRequest.RetryStrategy.SAFE_TO_RETRY), JsonNode.class);
    }

    // returns map of files for sample. The key is file ID be default, might be changes to names instead
    public String sampleLs(Integer sampleid, boolean byNames) {
        List filesList = DXSearch.findDataObjects().inFolder(DXContainer.getInstance(env.getProperty("dx-project")), "/samples/" + sampleid.toString() + "/rawdata").execute().asList();
        Map<String, Object> responseMap = new HashMap<>();
        String samplelsOutputResponse = null;

        if(filesList.size() > 0) {
            Map<String, Object> filesMap = new HashMap<>();
            DXFile dxFile;
            for (Object file : filesList) {
                dxFile = (DXFile) file;
                filesMap.putIfAbsent("fileId", dxFile.getId());
                filesMap.putIfAbsent("fileName", dxFile.describe().getName());
                filesMap.putIfAbsent("tags", dxFile.describe().getTags());
                if (byNames) {
                    responseMap.putIfAbsent(dxFile.describe().getName(), new HashMap<>(filesMap));
                } else {
                    responseMap.putIfAbsent(dxFile.getId(), new HashMap<>(filesMap));
                }
                filesMap.clear();
            }
        } 
        try {
                samplelsOutputResponse = new ObjectMapper().writeValueAsString(responseMap);
            } catch (JsonProcessingException e) {
                e.printStackTrace();
            }
        return samplelsOutputResponse;
    }
}
