package pl.intelliseq.genetraps.api.dx.helpers;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.ArrayNode;
import com.fasterxml.jackson.databind.node.ObjectNode;
import com.fasterxml.jackson.dataformat.yaml.YAMLFactory;
import com.mashape.unirest.http.HttpResponse;
import com.mashape.unirest.http.Unirest;
import com.mashape.unirest.http.exceptions.UnirestException;
import org.springframework.beans.factory.annotation.Value;

import javax.annotation.PostConstruct;
import java.io.*;
import java.net.URL;
import java.net.URLEncoder;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class WDLParserManager {

    @Value("${fileTree-prefix}")
    String pathPrefix;

    @Value("${fileTree-path}")
    String commonPath;

    @Value("${fileTree-ref}")
    String pathSuffix;

    @Value("${fileGet-prefix}")
    String filePathPrefix;

    @Value("${fileGet-ref}")
    String filePathRef;

    private JsonNode data;

    public JsonNode getData(){
        return data;
    }

    public JsonNode getData(String name){
        if (data.has(name)) {
            ObjectNode requestedData = new ObjectMapper().createObjectNode();
            requestedData.set(name, data.get(name));
            return requestedData;
        }
        else{
            return null;
        }
    }

    private JsonNode parserWDLToJSON(String path){
        String resourcePath = path;
        StringBuilder trueYaml = new StringBuilder();

        // Parsing file to clean YAML
        try {
            BufferedReader reader = new BufferedReader(new InputStreamReader(new URL(resourcePath).openStream()));

            // Omitting header
            for (int i=0;i<4;i++){
                reader.readLine();
            }

            // Parsing description to YAMLFactory
            String line;
            while ((line = reader.readLine()).startsWith("#")){
                if (line.length()>1){
                    line = line.substring(3);
                }
                else{
                    line = line.substring(1);
                }

                // if line is decoration
                if (line.length() >= 4 && (line.substring(2, 4).equals("--") || line.charAt(1) == '#')) {
                    continue;
                }
                trueYaml.append(line.replace("\\s","\\\\s")).append('\n');
            }

            String variableLine;
            Pattern descriptionPattern = Pattern.compile("#\\s*@Input.*");
            Pattern soloDescriptionPattern = Pattern.compile("((\\w*)\\s?=\\s?)((\\\".*?\\\")|(\\[.*\\])|(\\d*\\.?\\d?))");
            Pattern variablePatternQM = Pattern.compile(".*\\?.*");
            Pattern variablePatternEq = Pattern.compile("(.*?)\\s(.*)\\s*=\\s*(.*)");

            trueYaml.append("Inputs:\n");
            while((line = reader.readLine()) != null){

                // looking for /# @Input()/ lines
                if (descriptionPattern.matcher(line).matches()) {

                    variableLine = reader.readLine().trim();

                    // cutting variable line and parsing it to YAML
                    if (variablePatternQM.matcher(variableLine).matches()){
                        String[] parts = variableLine.split("\\? *");
                        trueYaml.append("  ").append(parts[1]).append(":\n");
                        trueYaml.append("    Type: ").append(parts[0]).append("\n");
                        trueYaml.append("    Optional: true\n");
                    }
                    else if (variablePatternEq.matcher(variableLine).matches()){
                        Matcher eq = variablePatternEq.matcher(variableLine);
                        eq.find();
                        trueYaml.append("  ").append(eq.group(2)).append(":\n");
                        trueYaml.append("    Type: ").append(eq.group(1)).append("\n");
                        trueYaml.append("    Default_Value: ").append(eq.group(3)).append("\n");
                        trueYaml.append("    Optional: false\n");
                    }
                    else
                    {
                        String[] parts = variableLine.split("\\s+");
                        trueYaml.append("  ").append(parts[1]).append(":\n");
                        trueYaml.append("    Type: ").append(parts[0]).append("\n");
                        trueYaml.append("    Optional: false\n");
                    }

                    // cutting found description line and parsing it to YAML
                    List<String> allMatches = new ArrayList<String>();
                    Matcher descriptionMatcher = soloDescriptionPattern.matcher(line);
                    while (descriptionMatcher.find())
                    {
                        allMatches.add(descriptionMatcher.group(2));
                        int i=6;
                        while(descriptionMatcher.group(i) == null) i--;
                        allMatches.add(descriptionMatcher.group(i));
                    }
                    for (int i=0;i<allMatches.size();i++)
                    {
                        trueYaml.append("    ").append(allMatches.get(i)).append(": ");
                        i++;

                        // if argument is list
                        if (allMatches.get(i).charAt(0) == '[')
                        {
                            String[] nodes = allMatches.get(i).substring(1, allMatches.get(i).length()-1).split(", ");
                            for (String node : nodes) {
                                trueYaml.append("\n    - ").append(node);
                            }
                            trueYaml.append('\n');
                        }
                        else {
                            trueYaml.append(allMatches.get(i)).append('\n');
                        }
                    }
                }
            }

            reader.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        //Parsing YAML to JSON
        ObjectMapper yamlReader = new ObjectMapper(new YAMLFactory());
        Object yamlObject;
        try {
            yamlObject = yamlReader.readValue(trueYaml.toString(), Object.class);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        ObjectMapper jsonWriter = new ObjectMapper();

        return jsonWriter.convertValue(yamlObject, JsonNode.class);
    }

    @PostConstruct
    public void collectData(){
        Pattern isWdlFile = Pattern.compile(".*\\.wdl");
        StringBuilder path;
        ObjectNode node = new ObjectMapper().createObjectNode();

        ArrayNode wdlDirectories = getListOfFiles((new StringBuilder(pathPrefix)).append(commonPath).append(pathSuffix).toString());

        // iterating through directories
        for (JsonNode fileNode: wdlDirectories)
        {
            String filePath = remover(fileNode.get("path").toString());
            path = (new StringBuilder(pathPrefix)).append(filePath).append(pathSuffix);

            ArrayNode wdlVersions = getListOfFiles(path.toString());

            // finding latest version in directory
            double highestVersion=0.0;
            for (JsonNode versionNode: wdlVersions){
                String versionPath = remover(versionNode.get("path").toString());
                double versionNumber = Double.parseDouble(remover(versionNode.get("name").toString()).substring(1));
                if (versionNumber>highestVersion) {
                    filePath = versionPath;
                    highestVersion=versionNumber;
                }
            }

            path = (new StringBuilder(pathPrefix)).append(filePath).append(pathSuffix);

            ArrayNode filesInVersionDirectory = getListOfFiles(path.toString());

            // finding and adding wdl file to node
            for (JsonNode file: filesInVersionDirectory){
                if (isWdlFile.matcher(remover(file.get("name").toString())).matches()){
                    String wdlPath;
                    try {
                        wdlPath = URLEncoder.encode(remover(file.get("path").toString()), StandardCharsets.UTF_8.toString());
                    } catch (UnsupportedEncodingException e) {
                        throw new RuntimeException(e);
                    }
                    StringBuilder apiPath = (new StringBuilder(filePathPrefix))
                            .append(wdlPath)
                            .append(filePathRef);
                    node.set(remover(file.get("name").toString()).substring(0, remover(file.get("name").toString()).length()-4), parserWDLToJSON(apiPath.toString()));
                }
            }
        }
        data = node;
    }

    private String remover(String oldString){
        return oldString.substring(1, oldString.length()-1);
    }

    private ArrayNode getListOfFiles(String path) {
        // Downloading JSON with WDL files data
        HttpResponse<String> response;
        try {
            response =  Unirest.get(path).asString();
        } catch (UnirestException e) {
            throw new RuntimeException(e);
        }

        // Parsing to ArrayNode
        ObjectMapper stackJson = new ObjectMapper();
        ArrayNode arrayOfFiles;
        try {
            arrayOfFiles = stackJson.readValue(response.getBody(), ArrayNode.class);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        return arrayOfFiles;
    }
}
