package pl.intelliseq.genetraps.api.dx.helpers;

import com.dnanexus.*;
import com.dnanexus.exceptions.ResourceNotFoundException;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.core.env.Environment;
import org.springframework.web.multipart.MultipartFile;
import pl.intelliseq.genetraps.api.dx.exceptions.PropertiesException;
import pl.intelliseq.genetraps.api.dx.exceptions.TagsException;
import pl.intelliseq.genetraps.api.dx.exceptions.WrongNumberOfFilesException;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.*;

public class DxApiProcessManager {

    @Autowired
    Environment env;

    public DXJob runUrlFetch(String inputUrl, Integer sampleId, String... tags) {
        var input = new HashMap<>();

        boolean leftTag = Arrays.asList(tags).contains("left");
        boolean rightTag = Arrays.asList(tags).contains("right");
        if(leftTag && rightTag) {
            throw new TagsException("both 'left' and 'right' tags cannot be used at same time");
        }
        if(!leftTag && !rightTag) {
            throw new TagsException("neither 'left' nor 'right' tag was used");
        }

        input.put("url", inputUrl);
        input.put("tags", tags);

        return DXApp.getInstance(env.getProperty("dx-url-fetch-ap-id"))
                .newRun()
                .setProject(DXProject.getInstance(env.getProperty("dx-project")))
                .setFolder(String.format("/samples/%s/rawdata", sampleId))
                .setInput(input)
                .run();
    }

    public String runUploadFile(MultipartFile mfile, Integer sampleId, String newfilename, List<String> tags) throws IOException {


        if(newfilename == null) {
            newfilename = mfile.getOriginalFilename();
        }

        DXFile.Builder filebuilder;
        filebuilder = DXFile.newFile()
                .setName(newfilename).setFolder(String.format("/samples/%s/rawdata", sampleId.toString()))
                .setProject(DXContainer.getInstance(env.getProperty("dx-project")));

        if(tags != null)
                filebuilder.addTags(tags);

        DXFile file;
        try {
            file = filebuilder.upload(mfile.getBytes()).build();
        }
        catch (ResourceNotFoundException e) {
            DXContainer.getInstance(env.getProperty("dx-project")).newFolder(String.format("/samples/%s/rawdata", sampleId.toString()), true);
            file = filebuilder.upload(mfile.getBytes()).build();
        }
        file.close();
        return file.getId();
    }

    public void runMkDir(Integer sampleId) {
        runMkDir(sampleId.toString());
    }

    public void runMkDir(String sampleId) {
        DXContainer.getInstance(env.getProperty("dx-project")).newFolder("/samples/" + sampleId, true);
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

    public DXJob runBwa(Integer sampleId) {
        var input = new HashMap<>();
        List fileListLeft = DXSearch.findDataObjects().inFolder(DXContainer.getInstance(env.getProperty("dx-project")), "/samples/" + sampleId.toString() + "/rawdata").withTag("left").execute().asList();
        if(fileListLeft.size() != 1) {
            throw new WrongNumberOfFilesException("expected: 1 with 'left' tag ; found: " + fileListLeft.size());
        }
        DXFile fastqFile1 = (DXFile) fileListLeft.get(0);

        List fileListRight = DXSearch.findDataObjects().inFolder(DXContainer.getInstance(env.getProperty("dx-project")), "/samples/" + sampleId.toString() + "/rawdata").withTag("right").execute().asList();
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

    public DXJob runGatkHC(Integer sampleId, String interval) {
        var input = new HashMap<>();
        DXFile bam = (DXFile) DXSearch.findDataObjects().nameMatchesRegexp(".*\\.bam$")
                .inFolder(DXContainer.getInstance(env.getProperty("dx-project")), "/samples/" + sampleId.toString() + "/bwa")
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

    public JsonNode JSONDescribe(String sampleId) {
        return DXJSON.safeTreeToValue(
                new DXHTTPRequest(DXEnvironment.create()).request("/" + sampleId + "/" + "describe",
                        new ObjectMapper().createObjectNode(), DXHTTPRequest.RetryStrategy.SAFE_TO_RETRY), JsonNode.class);
    }

    // returns map of files for sample. The key is file ID by default, might be changed to names instead
    public JsonNode sampleLs(Integer sampleId, boolean byNames) {
        List filesList = DXSearch.findDataObjects().inFolder(DXContainer.getInstance(env.getProperty("dx-project")), "/samples/" + sampleId.toString() + "/rawdata").execute().asList();
        Map<String, Object> responseMap = new HashMap<>();

        if(filesList.size() > 0) {
            HashMap<String, Object> filesMap = new HashMap<>();
            DXFile dxFile;
            for (Object file : filesList) {
                dxFile = (DXFile) file;
                filesMap.put("fileId", dxFile.getId());
                filesMap.put("fileName", dxFile.describe().getName());
                filesMap.put("tags", dxFile.describe().getTags());
                if (byNames) {
                    responseMap.put(dxFile.describe().getName(), new HashMap<>(filesMap));
                } else {
                    responseMap.put(dxFile.getId(), new HashMap<>(filesMap));
                }
                filesMap.clear();
            }
        }
        return new ObjectMapper().valueToTree(responseMap);
    }

    // creates and/or adds properties a specified sample folder
    public JsonNode propertiesPost(Integer sampleId, LinkedHashMap<String, String> properties) throws PropertiesException {
        List propertiesFileSearch = DXSearch.findDataObjects().nameMatchesExactly("properties").inFolder(DXContainer.getInstance(env.getProperty("dx-project")), String.format("/samples/%s", sampleId.toString())).execute().asList();
        DXFile file;
        if(propertiesFileSearch.size() == 0) {
            DXFile.Builder builder = DXFile.newFile().setName("properties")
                    .setFolder(String.format("/samples/%s", sampleId.toString()))
                    .putAllProperties(properties)
                    .setProject(DXContainer.getInstance(env.getProperty("dx-project")));
            try {
                file = builder.upload(new ByteArrayInputStream("properties".getBytes(StandardCharsets.UTF_8))).build().close();
            } catch(ResourceNotFoundException e) {
                throw new PropertiesException(String.format("cannot exist, because sample: %d does not exist", sampleId));
            }
        }
        else {
            file = (DXFile) propertiesFileSearch.get(0);
            Map<String, String> propertiesMap = file.describe(DXDataObject.DescribeOptions.get().withProperties()).getProperties();
            List<String> existingPropertiesMap = new LinkedList<>();
            int existingPropertiesCount = 0;
            for (String key: properties.keySet()) {
                if(propertiesMap.containsKey(key)) {
                    existingPropertiesCount++;
                    existingPropertiesMap.add(key);
                }
            }
            if(existingPropertiesCount != 0) {
                throw new PropertiesException(String.format("already contains properties of keys: %s", existingPropertiesMap.toString()));
            } else {
                file.putAllProperties(properties);
            }
        }
        return new ObjectMapper().valueToTree(file.describe(DXDataObject.DescribeOptions.get().withProperties()).getProperties());
    }

    // returns properties of a specified sample folder
    public JsonNode propertiesGet(Integer sampleId) throws PropertiesException {
        List propertiesFileSearch = DXSearch.findDataObjects().nameMatchesExactly("properties").inFolder(DXContainer.getInstance(env.getProperty("dx-project")), String.format("/samples/%s", sampleId.toString())).execute().asList();
        try {
            // dummy variable to tell if the folder exists -> throws an exception if it doesn't
            DXContainer.FolderContents x = DXContainer.getInstance(env.getProperty("dx-project")).listFolder(String.format("/samples/%s", sampleId.toString()));
        } catch(ResourceNotFoundException e) {
            throw new PropertiesException(String.format("cannot exist, because sample: %d does not exist", sampleId));
        }
        if(propertiesFileSearch.size() == 0) {
            throw new PropertiesException("does not exist");
        }
        DXFile file = (DXFile) propertiesFileSearch.get(0);
        return new ObjectMapper().valueToTree(file.describe(DXDataObject.DescribeOptions.get().withProperties()).getProperties());
    }

    // changes already existing properties of a specified sample folder
    public JsonNode propertiesPut(Integer sampleId, LinkedHashMap<String, String> properties) throws PropertiesException {
        List propertiesFileSearch = DXSearch.findDataObjects().nameMatchesExactly("properties").inFolder(DXContainer.getInstance(env.getProperty("dx-project")), String.format("/samples/%s", sampleId.toString())).execute().asList();
        try {
            // dummy variable to tell if the folder exists -> throws an exception if it doesn't
            DXContainer.FolderContents dummy = DXContainer.getInstance(env.getProperty("dx-project")).listFolder(String.format("/samples/%s", sampleId.toString()));
        } catch(ResourceNotFoundException e) {
            throw new PropertiesException(String.format("cannot exist, because sample: %d does not exist", sampleId));
        }
        DXFile file;
        if(propertiesFileSearch.size() == 0) {
            throw new PropertiesException(String.format("does not exist in sample: %d", sampleId));
        }
        else {
            file = (DXFile) propertiesFileSearch.get(0);
            Map<String, String> propertiesMap = file.describe(DXDataObject.DescribeOptions.get().withProperties()).getProperties();
            // will remember which properties couldn't be updated
            List<String> nonexistingPropertiesList = new LinkedList<>();
            // will remember the number of such properties
            int nonexistingPropertiesCount = 0;
            for (String key : properties.keySet()) {
                if (!propertiesMap.containsKey(key)) {
                    nonexistingPropertiesCount++;
                    nonexistingPropertiesList.add(key);
                }
            }
            if (nonexistingPropertiesCount != 0) {
                throw new PropertiesException(String.format("does not contain properties of keys: %s", nonexistingPropertiesList.toString()));
            } else {
                file.putAllProperties(properties);
            }
        }
        return new ObjectMapper().valueToTree(file.describe(DXDataObject.DescribeOptions.get().withProperties()).getProperties());
    }

    // deletes already existing properties of a specified sample folder
    public JsonNode propertiesDelete(Integer sampleId, LinkedHashMap<String, String> properties) throws PropertiesException {
        List propertiesFileSearch = DXSearch.findDataObjects().nameMatchesExactly("properties").inFolder(DXContainer.getInstance(env.getProperty("dx-project")), String.format("/samples/%s", sampleId.toString())).execute().asList();
        try {
            // dummy variable to tell if the folder exists -> throws an exception if it doesn't
            DXContainer.FolderContents dummy = DXContainer.getInstance(env.getProperty("dx-project")).listFolder(String.format("/samples/%s", sampleId.toString()));
        } catch(ResourceNotFoundException e) {
            throw new PropertiesException(String.format("cannot exist, because sample: %d does not exist", sampleId));
        }
        DXFile file;
        if(propertiesFileSearch.size() == 0) {
            throw new PropertiesException(String.format("does not exist in sample: %d", sampleId));
        }
        else {
            file = (DXFile) propertiesFileSearch.get(0);
            Map<String, String> propertiesMap = file.describe(DXDataObject.DescribeOptions.get().withProperties()).getProperties();
            // will remember which properties couldn't be deleted
            List<String> nonexistingPropertiesList = new LinkedList<>();
            // will remember the number of such properties
            int nonexistingPropertiesCount = 0;
            for (String key : properties.keySet()) {
                if (!propertiesMap.containsKey(key)) {
                    nonexistingPropertiesCount++;
                    nonexistingPropertiesList.add(key);
                }
            }
            if (nonexistingPropertiesCount != 0) {
                throw new PropertiesException(String.format("does not contain properties of keys: %s", nonexistingPropertiesList.toString()));
            } else {
                for (String key : properties.keySet()) {
                    file.removeProperty(key);
                }
            }
        }
        return new ObjectMapper().valueToTree(file.describe(DXDataObject.DescribeOptions.get().withProperties()).getProperties());
    }
}
