package pl.intelliseq.genetraps.api.dx.parser;

public class DxJob {

	private String name;
	private String state;
	private String id;
	private String owner;
	private String date;
	private String runtime;
	
	public static DxJob getDxJobById(String jobId) {
		String result = DxRunner.runCommand("dx find jobs --id " + jobId);
		
		return new DxJob(result);
	}
	
	private DxJob(String dxRunnerResult) {
		String[] elements = dxRunnerResult.replace("\n", "").split(" ");
		this.name = elements[1];
		this.state = removeFirstAndLastChar(elements[2]);
		this.id = elements[3];
		this.owner = elements[5];
		this.date = elements[6] + " " + elements[7];
		if (elements.length > 9)
			this.runtime = removeLastChar(elements[9]);
		if (elements.length > 10)
			this.runtime = removeLastChar(elements[10]);
	}
	
	private static String removeFirstAndLastChar(String str) {
	    return str.substring(1, str.length() - 1);
	}
	
	private static String removeLastChar(String str) {
	    return str.substring(0, str.length() - 1);
	}
	
	
	@Override
	public String toString() {
		return name + " :: " + state + " :: " + id + " :: " + owner + " :: " + date + " :: " + runtime;
	}
	
}
