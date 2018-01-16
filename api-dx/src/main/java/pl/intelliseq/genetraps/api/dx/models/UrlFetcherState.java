package pl.intelliseq.genetraps.api.dx.models;


import pl.intelliseq.genetraps.api.dx.enums.Paths;

import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * Created by ltw on 07.06.17.
 */
public class UrlFetcherState {

    private String jobId;
    private String state;
    private String failureMessage;

    public UrlFetcherState(String jobId)  {

        this.jobId = jobId;

        List<String> urlFetcherState = GetUrlFetcherState(this.jobId);

        this.state = urlFetcherState.get(0);
        this.failureMessage = urlFetcherState.get(1);
    }

    public List<String> GetUrlFetcherState(String jobId) {

        String state = null;
        String failureMessage = null;
        String scriptPath = Paths.GET_URL_FETCHER_STATE.getPath();

        ProcessBuilder pb = new ProcessBuilder(scriptPath, jobId);
        Map<String, String> env = pb.environment();

        try {
            Process p = pb.start();

            String stdOut = getInputAsString(p.getInputStream());
            String stdErr = getInputAsString(p.getErrorStream());

            System.out.println(stdOut);
            System.out.println(stdErr);

            String stdOutLines[] = stdOut.split("\\r?\\n");



            try {
                state = stdOutLines[0].substring(stdOut.lastIndexOf("State: ") + 7).replaceAll("\\s+", "");
            } catch (StringIndexOutOfBoundsException e) {
                e.printStackTrace();
            }

            if ( stdOutLines.length == 2) {

                try {
                    failureMessage = stdOutLines[1].substring(stdOut.lastIndexOf("Failure message: ") + 17).replaceAll("\\s+", "");
                } catch (StringIndexOutOfBoundsException e) {
                    e.printStackTrace();
                }

            }

        } catch (IOException e) {
            e.printStackTrace();
        }


        return Arrays.asList(state, failureMessage);

    }

    public String getJobId() {
        return jobId;
    }
    public String getState() {
        return state;
    }
    public String getFailureMessage() { return failureMessage; }


    private String getInputAsString(InputStream is)
    {
        try(java.util.Scanner s = new java.util.Scanner(is))
        {
            return s.useDelimiter("\\A").hasNext() ? s.next() : "";
        }
    }


}
