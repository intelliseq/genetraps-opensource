import { Component, Input, Inject, Injectable } from '@angular/core';
import { Http, Headers } from '@angular/http';
import 'rxjs/operator/map';
import 'rxjs/Rx';
import { PhenotypeService } from './phenotype.service';
import { ActivatedRoute } from '@angular/router';

import { APP_CONFIG, IAppConfig } from './app.config';

@Component({
  selector: 'my-similar-diseases',
  template: `
  <div class="row pb-1">
    <div class="col-xs-12">
    <h4>Input disease</h4>
        <input  ng2-auto-complete
          [(ngModel)]="diseaseName"
          (ngModelChange)="getPhenotypes($event)"
          placeholder="Enter Disease Name"
          [source]="diseaseAutocompleteAdress"
          path-to-data="results"
          min-chars="0" />
    </div>
  </div>
  <div class="row pb-1">
    <div class="col-xs-12">
    <h4>{{diseaseName}}</h4>
    </div>
  </div>
  `,
  styles: [`
    ng2-auto-complete, input {
      display: block; border: 0px solid #ccc; width: 600px;
    }
  `],
})
export class SimilarDiseasesComponent {
  private diseaseName : string;

  getPhenotypes(disease : string) {
    this.phenotypeService.getPhenotypesByDisease(disease);
  }

  diseaseAutocompleteAdress: string = this.config.apiEndpoint + "/disease-autocomplete?firstLetters=:keyword&resultsCount=20";

  constructor(private phenotypeService: PhenotypeService,
  @Inject(APP_CONFIG) private config: IAppConfig) {}

}
