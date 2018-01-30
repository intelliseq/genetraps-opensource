package pl.intelliseq.genetraps.api.dx.models;

import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import pl.intelliseq.genetraps.api.dx.exceptions.IseqParseException;

public class IseqJSON {
    private JSONObject jsonObject;

    public IseqJSON(String string){
        JSONParser jsonParser = new JSONParser();
        try {
            jsonObject = (JSONObject) jsonParser.parse(string);
        } catch (ParseException e) {
            throw new IseqParseException(e);
        }
    }

    public IseqJSON(JSONObject jsonObject){
        this.jsonObject = jsonObject;
    }

    public IseqJSON(String key, String value){
        jsonObject = new JSONObject();
        jsonObject.put(key, value);
    }

    public IseqJSON put(Object key, Object value){
        jsonObject.put(key, value);
        return this;
    }

    public String getString(String key){
        return (String) jsonObject.get(key);
    }

    public IseqJSON getIseqJSON(String key){
        return new IseqJSON((JSONObject) jsonObject.get(key));
    }

    public Object get(String key){
        return jsonObject.get(key);
    }

    @Override
    public String toString() {
        return jsonObject.toJSONString();
    }
}
