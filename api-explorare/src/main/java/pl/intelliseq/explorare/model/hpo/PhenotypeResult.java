package pl.intelliseq.explorare.model.hpo;

import java.util.List;

import com.fasterxml.jackson.annotation.JsonView;

import pl.intelliseq.explorare.utils.json.Views;

public class PhenotypeResult {
	
	@JsonView(Views.Rest.class)
	private HpoTerm phenotype;
	
	@JsonView(Views.Rest.class)
	private HpoTerm parent;
	
	@JsonView(Views.Rest.class)
	private List<HpoTerm> children;
	
	public PhenotypeResult(HpoTerm phenotype) {
		this.phenotype = phenotype;
		this.parent = phenotype.getParent();
		this.children = phenotype.getChildren();
	}

	public HpoTerm getPhenotype() {
		return phenotype;
	}

	public void setPhenotype(HpoTerm phenotype) {
		this.phenotype = phenotype;
	}

	public HpoTerm getParent() {
		return parent;
	}

	public void setParent(HpoTerm parent) {
		this.parent = parent;
	}

	public List<HpoTerm> getChildren() {
		return children;
	}

	public void setChildren(List<HpoTerm> children) {
		this.children = children;
	}

}
