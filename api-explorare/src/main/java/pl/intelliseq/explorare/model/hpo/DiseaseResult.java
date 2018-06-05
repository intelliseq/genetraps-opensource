package pl.intelliseq.explorare.model.hpo;

public class DiseaseResult implements Comparable {
	private String name;
	private Double score;
	
	public DiseaseResult(String name, Double score) {
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
		int compareResult = -this.score.compareTo(((DiseaseResult) o).getScore());
		if (compareResult != 0) return compareResult;
		else return this.name.compareTo(((DiseaseResult) o).getName());
	}
}
