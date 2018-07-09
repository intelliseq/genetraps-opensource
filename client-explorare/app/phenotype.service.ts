import { Injectable, Inject, OnInit } from '@angular/core';
import { BehaviorSubject }    from 'rxjs/BehaviorSubject';
import { Http, Headers } from '@angular/http';
import 'rxjs/operator/map';
import 'rxjs/Rx';
import { Phenotype } from './phenotype';
import { Disease } from './disease';
import { Gene } from './gene';

import { APP_CONFIG, IAppConfig } from './app.config';

@Injectable()
export class PhenotypeService {

  /*** parsed text ***/
  parsedText = new BehaviorSubject<string>("");

  /*** phenotypes ***/
  phenotypes = new BehaviorSubject<Phenotype[]>([]);
  diseases = new BehaviorSubject<Disease[]>([]);
  genes = new BehaviorSubject<Gene[]>([]);

  /*** getters & setters ***/
  addPhenotype(phenotype : Phenotype) {
    var updatedPhenotypes : Phenotype[] = this.phenotypes.value;
    updatedPhenotypes.push(phenotype);
    this.phenotypes.next(updatedPhenotypes);
    //this.phenotypes;
  }

  /*** requests ***/
  parseText(query: string) {
    var data = JSON.stringify({ "query": query });
    var headers = new Headers();
    headers.append('Content-Type', 'application/json');
    //console.log(this.query);
    this.http.post(this.config.apiEndpoint + '/parse-text', data, {headers: headers})
              .map(response => response.json())
              .subscribe(result =>
                { this.parsedText.next(result.markedText);
                  this.phenotypes.next(<Phenotype[]> result.hpoTerms);
                  console.log(result.hpoTerms);});
  }

  getPhenotypesByDisease(disease: string) {
    var parsedDisease = disease.replace(/ /g, '%20');
    var url = this.config.apiEndpoint + "/get-phenotypes-by-disease?disease="+parsedDisease;
    this.http.get(url)
              .map(response => response.json())
              .subscribe(result =>
                { this.phenotypes.next(<Phenotype[]> result.results[0]);
                  console.log(result.hpoTerms);});
  }

  getDiseasesByPhenotypes(phenotypes : Phenotype[]) {
    //console.log("disease request");
    var data = "{ \"hpoTerms\": " + JSON.stringify(phenotypes) + " }";
    var headers = new Headers();
    headers.append('Content-Type', 'application/json');
    //console.log(this.query);
    this.http.post(this.config.apiEndpoint + '/get-diseases', data, {headers: headers})
              .map(response => response.json())
              .subscribe(result =>
                { //this.parsedText.next(result.markedText);
                  this.diseases.next(<Disease[]> result);
                  //console.log(result);
                  console.log("disease request completed");});
  }

  getGenesByPhenotypes(phenotypes : Phenotype[]) {
    //console.log("disease request");
    var data = "{ \"hpoTerms\": " + JSON.stringify(phenotypes) + " }";
    var headers = new Headers();
    headers.append('Content-Type', 'application/json');
    //console.log(this.query);
    this.http.post(this.config.apiEndpoint + '/get-genes', data, {headers: headers})
              .map(response => response.json())
              .subscribe(result =>
                { //this.parsedText.next(result.markedText);
                  this.genes.next(<Disease[]> result);
                  //console.log(result);
                  console.log("gene request completed");});
  }

  constructor( private http: Http,
  @Inject(APP_CONFIG) private config: IAppConfig ) {
    this.phenotypes.subscribe(value => this.getDiseasesByPhenotypes(value));
    this.phenotypes.subscribe(value => this.getGenesByPhenotypes(value));
  }

}
