package pl.intelliseq.genetraps.api.dx.helpers;

import com.dnanexus.*;
import com.dnanexus.exceptions.WrongNumberOfFilesException;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.core.env.Environment;

import java.util.HashMap;
import java.util.List;

public class DxApiProcessManager {

    @Autowired
    Environment env;

    public DXJob runUrlFetch(String inputUrl, String sampleNumber, String... tags) {
        var input = new HashMap<>();
        input.put("url", inputUrl);
        input.put("tags", tags);

        return DXApp.getInstance(env.getProperty("dx-url-fetch-ap-id"))
                .newRun()
                .setProject(DXProject.getInstance(env.getProperty("dx-project")))
                .setFolder(String.format("/samples/%s/rawdata", sampleNumber))
                .setInput(input)
                .run();
    }

    public void runMkDir(Integer sampleNumber) {
        runMkDir(sampleNumber.toString());
    }

    public void runMkDir(String sampleNumber) {
        DXContainer.getInstance(env.getProperty("dx-project")).newFolder("/samples/" + sampleNumber, true);
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

    public DXJob runBwa(String fastq_fileId_1, String fastq_fileId_2, String reference) {
        var input = new HashMap<>();
        DXFile fastqFile1 = DXFile.getInstance(fastq_fileId_1);
        DXFile fastqFile2 = DXFile.getInstance(fastq_fileId_2);
        DXFile ref = DXFile.getInstance(reference);
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

    public DXJob runBwa(Integer samples_number, String reference) {
        var input = new HashMap<>();
        List fileListLeft = DXSearch.findDataObjects().inFolder(DXContainer.getInstance(env.getProperty("dx-project")), "/samples/" + samples_number.toString() + "/rawdata").withTag("left").execute().asList();
        if(fileListLeft.size() != 1) {
            throw new WrongNumberOfFilesException("expected: 1 with 'left' tag ; found: " + fileListLeft.size());
        }
        DXFile fastqFile1 = (DXFile) fileListLeft.get(0);

        List fileListRight = DXSearch.findDataObjects().inFolder(DXContainer.getInstance(env.getProperty("dx-project")), "/samples/" + samples_number.toString() + "/rawdata").withTag("right").execute().asList();
        if(fileListRight.size() != 1) {
            throw new WrongNumberOfFilesException("expected: 1 with 'right' tag ; found: " + fileListRight.size());
        }
        DXFile fastqFile2 = (DXFile) fileListRight.get(0);

        DXFile ref = DXFile.getInstance(reference);
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

    public DXApplet getAppletFromName(String appletName) {
        return (DXApplet) DXSearch.findDataObjects().nameMatchesExactly(appletName).execute().asList().get(0);
    }

    public JsonNode JSONDescribe(String id) {
        return DXJSON.safeTreeToValue(
                new DXHTTPRequest(DXEnvironment.create()).request("/" + id + "/" + "describe",
                        new ObjectMapper().createObjectNode(), DXHTTPRequest.RetryStrategy.SAFE_TO_RETRY), JsonNode.class);
    }

}
