package pl.intelliseq.genetraps.api.dx.parser;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;

public class DxRunner {

	private static Logger log = Logger.getLogger(DxRunner.class);
	
	/*
	 * example: "printf \"Y\\n\" | dx run touch -iname=test"
	 */
	public static String runCommand(String command) {
		
		log.info("Running: " + command);
		
		ProcessBuilder builder = new ProcessBuilder();
    	builder.command("sh", "-c", command);
    	builder.directory(new File(System.getProperty("user.home")));
    	
    	log.info("Builder: " + builder);
    	
    	Process process;
		try {
			process = builder.start();
			log.info("Process: " + process);
		} catch (IOException e) {
			throw new DxRunnerException(e);
		}
    	String result = new BufferedReader(new InputStreamReader(process.getInputStream()))
    			  .lines().collect(Collectors.joining("\n"));
    	log.info("Result: " + result);
    	
    	return result;
	}
	
}
