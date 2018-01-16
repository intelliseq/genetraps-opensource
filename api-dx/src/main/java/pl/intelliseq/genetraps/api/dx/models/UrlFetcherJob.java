package pl.intelliseq.genetraps.api.dx.models;

import pl.intelliseq.genetraps.api.dx.enums.Paths;

import java.io.IOException;
import java.io.InputStream;
import java.util.Map;

/**
 * Created by ltw on 05.06.17.
 */
public class UrlFetcherJob {

    private String url;
    private String id;

    public UrlFetcherJob(String url)  {
        this.url = url;
        this.id = RunUrlFetcherApp(this.url);
    }

    public String RunUrlFetcherApp(String url) {

        String id = null;
        String scriptPath = Paths.RUN_URL_FETCHER.getPath();

        ProcessBuilder pb = new ProcessBuilder(scriptPath, url);
        Map<String, String> env = pb.environment();
        
        try {
            Process p = pb.start();

            String stdOut = getInputAsString(p.getInputStream());
            String stdErr = getInputAsString(p.getErrorStream());

            System.out.println(stdOut);
            System.out.println(stdErr);

            try {
                id = stdOut.substring(stdOut.lastIndexOf("Job ID: ") + 8).replaceAll("\\s+", "");
            } catch (StringIndexOutOfBoundsException e) {

                System.out.println(stdOut);
                e.printStackTrace();
            }

        } catch (IOException e) {
            e.printStackTrace();
        }

        return id;

    }

    public String getUrl() {
        return url;
    }
    public String getId() {
        return id;
    }


    private String getInputAsString(InputStream is)
    {
        try(java.util.Scanner s = new java.util.Scanner(is))
        {
            return s.useDelimiter("\\A").hasNext() ? s.next() : "";
        }
    }



}
