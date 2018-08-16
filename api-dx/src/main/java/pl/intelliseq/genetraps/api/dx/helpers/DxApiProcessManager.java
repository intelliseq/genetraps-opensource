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
        input.put("fastq_file1", fastqFile1);
        input.put("fastq_file2", fastqFile2);
        input.put("reference", reference);


        return getAppletFromName("iseq_bwa")
                .newRun()
                .setProject(DXProject.getInstance(env.getProperty("dx-project")))
                .setInput(input)
                .setFolder(fastqFile1.describe().getFolder().replace("rawdata", "bwa/" + fastqFile1.describe().getName()))
                .setFolder(fastqFile2.describe().getFolder().replace("rawdata", "bwa/" + fastqFile2.describe().getName()))
                .run();
    }

    public DXJob runBwa(Integer samples_number, String reference) {
        var input = new HashMap<>();
        List fileList = DXSearch.findDataObjects().inFolder(DXContainer.getInstance(env.getProperty("dx-project")), "/samples/" + samples_number.toString()).withTag("left").execute().asList();
        if(fileList.size() != 1) {
            throw new WrongNumberOfFilesException("expected: 1 with 'left' tag ; found: " + fileList.size());
        }
        DXFile fastqFile1 = (DXFile) fileList.get(0);
        fileList = DXSearch.findDataObjects().withTag("right").execute().asList();
        if(fileList.size() != 1) {
            throw new WrongNumberOfFilesException("expected: 1 with 'right' tag ; found: " + fileList.size());
        }
        DXFile fastqFile2 = (DXFile) fileList.get(0);
        input.put("fastq_file1", fastqFile1);
        input.put("fastq_file2", fastqFile2);
        input.put("reference", reference);


        return getAppletFromName("iseq_bwa")
                .newRun()
                .setProject(DXProject.getInstance(env.getProperty("dx-project")))
                .setInput(input)
                .setFolder(fastqFile1.describe().getFolder().replace("rawdata", "bwa/" + fastqFile1.describe().getName()))
                .setFolder(fastqFile2.describe().getFolder().replace("rawdata", "bwa/" + fastqFile2.describe().getName()))
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
