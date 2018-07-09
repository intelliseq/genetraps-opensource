import { Component, Input, Injectable } from '@angular/core';
import { Http, Headers } from '@angular/http';
import 'rxjs/operator/map';
import 'rxjs/Rx';
import { PhenotypeService } from './phenotype.service';
import { ActivatedRoute } from '@angular/router';

@Component({
  selector: 'my-pheno-marks',
  template: `
  <div class="row pb-1">
    <div class="col-sm-9">
      <textarea [(ngModel)]="query" rows="5" class="border-no width-full" placeholder="Paste epicrysis here" (change)=postItems()></textarea>
    </div>
    <div class="col-sm-3">
      <a class="btn btn-info width-full mb-1" style="color:white;" (click)=example()>example</a>
      <a class="btn btn-success width-full mb-1" style="color:white;" (click)=postItems()>submit</a>
    </div>
  </div>
  <div class="row pb-1">
    <div [innerHTML]="output" style="clear: both;"></div>
  </div>

  `,
})
export class PhenoMarksComponent {

  private EXAMPLE_MICRO: string = "Warburg et al. (1993) used the designation Micro syndrome for an autosomal recessive syndrome comprising microcephaly, microcornea, congenital cataract, mental retardation, optic atrophy, and hypogenitalism. They described an affected brother and sister and their male cousin. The sibs were offspring of a consanguineous Pakistani marriage; the parents of the cousin denied consanguinity. Agenesis of the corpus callosum, prominent root of the nose, large anteverted ears, facial hypertrichosis, small pupils with posterior synechiae, hypotonia, mild to moderate spastic palsy with hip dislocations, and hormonal dysfunction, presumably of hypothalamic origin, were other features. The children were almost blind, whether or not the cataracts had been operated on. The electroretinographic responses indicated dysfunction of both retinal rods and cones, and the visual evoked potentials confirmed optic nerve atrophy. The children were late walkers and were incontinent of urine and stools. In the differential diagnosis, Warburg et al. (1993) considered COFS syndrome (214150), CAMAK/CAMFAK syndromes (212540), Martsolf syndrome (212720), lethal Rutledge syndrome (270400), and lethal Neu-Laxova syndrome (256520).";
  private query: string = "";
  private output: string = "";

  example() {
    this.query = this.EXAMPLE_MICRO;
  }

  postItems() {
    this.output = `
      <div class="spinner">
        <div class="bounce1"></div>
        <div class="bounce2"></div>
        <div class="bounce3"></div>
      </div>`;
    this.phenotypeService.parseText(this.query);
  }

  constructor(private phenotypeService: PhenotypeService) {
    phenotypeService.parsedText.subscribe(value => this.output = value);
  }

}
