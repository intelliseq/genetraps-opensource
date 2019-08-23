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
        if (leftTag && rightTag) {
            throw new TagsException("both 'left' and 'right' tags cannot be used at same time");
        }
        if (!leftTag && !rightTag) {
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

    // creates and/or adds properties a specified sample folder
    public JsonNode propertiesPost(Integer sampleId, LinkedHashMap<String, String> properties) throws PropertiesException {
        List propertiesFileSearch = DXSearch.findDataObjects().nameMatchesExactly("properties").inFolder(DXContainer.getInstance(env.getProperty("dx-project")), String.format("/samples/%s", sampleId.toString())).execute().asList();
        DXFile file;
        if (propertiesFileSearch.size() == 0) {
            DXFile.Builder builder = DXFile.newFile().setName("properties")
                    .setFolder(String.format("/samples/%s", sampleId.toString()))
                    .putAllProperties(properties)
                    .setProject(DXContainer.getInstance(env.getProperty("dx-project")));
            try {
                file = builder.upload(new ByteArrayInputStream("properties".getBytes(StandardCharsets.UTF_8))).build().close();
            } catch (ResourceNotFoundException e) {
                throw new PropertiesException(String.format("cannot exist, because sample: %d does not exist", sampleId));
            }
        } else {
            file = (DXFile) propertiesFileSearch.get(0);
            Map<String, String> propertiesMap = file.describe(DXDataObject.DescribeOptions.get().withProperties()).getProperties();
            List<String> existingPropertiesMap = new LinkedList<>();
            int existingPropertiesCount = 0;
            for (String key : properties.keySet()) {
                if (propertiesMap.containsKey(key)) {
                    existingPropertiesCount++;
                    existingPropertiesMap.add(key);
                }
            }
            if (existingPropertiesCount != 0) {
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
        } catch (ResourceNotFoundException e) {
            throw new PropertiesException(String.format("cannot exist, because sample: %d does not exist", sampleId));
        }
        if (propertiesFileSearch.size() == 0) {
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
        } catch (ResourceNotFoundException e) {
            throw new PropertiesException(String.format("cannot exist, because sample: %d does not exist", sampleId));
        }
        DXFile file;
        if (propertiesFileSearch.size() == 0) {
            throw new PropertiesException(String.format("does not exist in sample: %d", sampleId));
        } else {
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
        } catch (ResourceNotFoundException e) {
            throw new PropertiesException(String.format("cannot exist, because sample: %d does not exist", sampleId));
        }
        DXFile file;
        if (propertiesFileSearch.size() == 0) {
            throw new PropertiesException(String.format("does not exist in sample: %d", sampleId));
        } else {
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
