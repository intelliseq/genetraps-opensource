package pl.intelliseq.genetraps.api.dx.helpers;

import org.apache.log4j.Logger;
import pl.intelliseq.genetraps.api.dx.exceptions.DxRunnerException;
import pl.intelliseq.genetraps.api.dx.models.IseqJSON;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.time.Duration;
import java.time.Instant;
import java.util.stream.Collectors;

public class ProcessManager {
    private Logger log = Logger.getLogger(ProcessManager.class);

    private void dumpProcessInfo(ProcessHandle ph){
        log.debug("PROCESS INFORMATION");
        log.debug("===================");
        log.debug(String.format("Process id: %d%n", ph.pid()));
        ProcessHandle.Info info = ph.info();
        log.debug(String.format("Command: %s%n", info.command().orElse("")));
        String[] args = info.arguments().orElse(new String[]{});
        log.debug("Arguments:");
        for (String arg: args)
            log.debug(String.format("   %s%n", arg));
        log.debug(String.format("Command line: %s%n", info.commandLine().orElse("")));
        log.debug(String.format("Start time: %s%n",
                info.startInstant().orElse(Instant.now()).toString()));
        log.debug(String.format("Run time duration: %sms%n",
                info.totalCpuDuration()
                        .orElse(Duration.ofMillis(0)).toMillis()));
        log.debug(String.format("Owner: %s%n", info.user().orElse("")));
    }

    public String runCommand(String command) {
        ProcessBuilder builder = new ProcessBuilder();
        log.info("Command:\t"+command);
        builder.command("sh", "-c", command);
        Process p;
        String result;
        try {
            p = builder.start();
//            dumpProcessInfo(p.toHandle());
            result = new BufferedReader(new InputStreamReader(p.getInputStream()))
                    .lines().collect(Collectors.joining("\n"));
//            dumpProcessInfo(p.toHandle());
        } catch (IOException e) {
            throw new DxRunnerException(e);
        }
        return result;
    }

    public String runCommandAndGetJobId(String command){
        return runCommand(command+" -y | tail -1 | grep -Eo \"job-\\w{24}\"");
    }

    public String runTouch(String args){
        return runCommand("dx run touch "+args);
    }

    public String runMkdir(String args){
        return runCommand("dx mkdir samples/"+args);
    }

    public IseqJSON runJSONDescribe(String args){
        return new IseqJSON(runCommand("dx describe --json "+args));
    }

    public String runFastqc(String args){
        String command = String.format("dx run iseq-fastqc -i fastq_file=\"%s\"", args);
        return runCommandAndGetJobId(command);
    }

    public String runUrlFetch(String inputUrl, String sampleNumber, String... tags){
        //TODO: tag
        log.info(tags);
        StringBuilder command = new StringBuilder("dx run url_fetcher");
        for(String tag:tags){
            command.append(" -i tags=").append(tag);
        }
        command.append(String.format(" -i url=\"%s\" --folder=\"samples/%s/rawdata\"", inputUrl, sampleNumber));
        return runCommandAndGetJobId(command.toString());
    }

    public String runBwa(String left, String right, String outputFolder){
        String command = String.format("dx run bwa_mem_fastq_read_mapper -i reads_fastqgz=%s -i reads2_fastqgz=%s -i genomeindex_targz=project-BQpp3Y804Y0xbyG4GJPQ01xv:file-BFBy4G805pXZKqV1ZVGQ0FG8 --destination %s", left, right, outputFolder);
        return runCommandAndGetJobId(command);
    }
}
