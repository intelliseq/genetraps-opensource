package pl.intelliseq.explorare;

import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;

import pl.intelliseq.explorare.model.hpo.DiseaseGeneDictionary;
import pl.intelliseq.explorare.model.hpo.HpoOboParser;
import pl.intelliseq.explorare.model.hpo.HpoTree;
import pl.intelliseq.explorare.model.phenoMarks.PhenoMarksParser;

@SpringBootApplication
@Configuration
public class ExplorareServer {

	HpoTree hpoTree = new HpoOboParser().getHpoTree();
	
    public static void main(String[] args) {
        SpringApplication.run(ExplorareServer.class, args);
    }

    @Bean
    PhenoMarksParser medicMarksParser() {
        return new PhenoMarksParser(this.hpoTree);
    }
    
    @Bean
    HpoTree hpoTree() {
    	return this.hpoTree;
    }
    
    @Bean
    DiseaseGeneDictionary diseaseGeneDictionary() {
    	return new DiseaseGeneDictionary();
    }
    
}