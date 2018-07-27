package pl.intelliseq.explorare.model.hpo;

public class GeneResult implements Comparable {
	private String name;
	private Double score;
	
	public GeneResult(String name, Double score) {
		this.name = name;
		this.score = score;
	}
	
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public Double getScore() {
		return score;
	}
	public void setScore(Double score) {
		this.score = score;
	}

	@Override
	public int compareTo(Object o) {
		if (this.name.equals(((GeneResult) o).getName())) {System.out.println("Compare " + this.name);return 0;}
		int compareResult = -this.score.compareTo(((GeneResult) o).getScore());
		if (compareResult != 0) return compareResult;
		else return this.name.compareTo(((GeneResult) o).getName());
	}
}
