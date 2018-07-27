package pl.intelliseq.explorare.model;

import java.util.HashSet;
import java.util.Set;

import com.fasterxml.jackson.annotation.JsonView;

import pl.intelliseq.explorare.utils.json.Views;

public class Results {
	
	@JsonView(Views.Rest.class)
	Set <Object> results = new HashSet <Object>();

	public Set<Object> getResults() {
		return results;
	}

	public void addResult(Object result) {
		this.results.add(result);
	}
	
	public boolean sizeGreaterOrEqualTo(int count) {
		return results.size() >= count;
	}
	
}
