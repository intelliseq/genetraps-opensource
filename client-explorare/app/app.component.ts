import { Component } from '@angular/core';
import { PhenotypeService } from './phenotype.service';

@Component({
  selector: 'my-app',
  template: `
<div class="container">
  <div class="row py-1" data-spy="affix">
    <div class="col-sm-4">
      <h1>explorare</h1>
    </div>
    <div class="col-sm-2">
      Phenotypes <span class="badge badge-pill badge-success">{{phenotypeCount}}</span>
    </div>
    <div class="col-sm-2">
      Diseases <span class="badge badge-pill badge-success">{{diseaseCount}}</span>
    </div>
    <div class="col-sm-2">
      Genes <span class="badge badge-pill badge-success">{{geneCount}}</span>
    </div>
  </div>
  <div class="row pb-1">

    <div class="col-sm-3">
        <a class="btn btn-primary btn-block px-0 mb-1" style="white-space: normal;" [routerLink]="['/pheno-marks',{dataToAdd:'a'}]" routerLinkActive="active">
        Pheno Marks</a>
        <p class="small text-center">paste epicrysis and get phenotypes</p>
    </div>
    <div class="col-sm-3">
        <a class="btn btn-primary btn-block px-0 mb-1" style="white-space: normal;" routerLink="/similar-diseases" routerLinkActive="active">
        Disease Input</a>
        <p class="small text-center">input disease and get phenotypes</p>
    </div>
    <div class="col-sm-3">
        <a class="btn btn-primary btn-block px-0 mb-1" style="white-space: normal;" routerLink="/pheno-panel" routerLinkActive="active">
        Pheno Panel</a>
        <p class="small text-center">curate or enter phenotypes</p>
    </div>
    <!--<div class="col-6 col-sm-2">
        <a class="btn btn-primary btn-block px-0 mb-1" style="white-space: normal;" routerLink="/pheno-panel" routerLinkActive="active">Genes</a>
    </div>
    <div class="col-6 col-sm-2">
        <a class="btn btn-primary btn-block px-0 mb-1" style="white-space: normal;" routerLink="/pheno-panel" routerLinkActive="active">Diseases</a>
    </div>-->
  </div>
  <router-outlet></router-outlet>
</div>
  `,
})
export class AppComponent  {
  phenotypeCount: number;
  diseaseCount: number;
  geneCount: number;
  constructor (phenotypeService: PhenotypeService) {
    phenotypeService.phenotypes.subscribe(value => {this.phenotypeCount = value.length;});
    phenotypeService.diseases.subscribe(value => {this.diseaseCount = value.length;});
    phenotypeService.genes.subscribe(value => {this.geneCount = value.length;});
  }
}
