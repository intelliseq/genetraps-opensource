import { Component, Input, Inject, Injectable, ChangeDetectorRef } from '@angular/core';
import { BehaviorSubject }    from 'rxjs/BehaviorSubject';
import { Http, Headers } from '@angular/http';
import 'rxjs/operator/map';
import 'rxjs/Rx';

import { PhenotypeService } from './phenotype.service';
import { Phenotype } from './phenotype';
import { Disease } from './disease';
import { Gene } from './gene';

import { ViewContainerRef } from '@angular/core';

import { APP_CONFIG, IAppConfig } from './app.config';

@Component({
  selector: 'my-pheno-panel',
  template: `


  <div class="row pb-1">
  <div class="col-sm-3">
      add phenotype<input  ng2-auto-complete
        [(ngModel)]="phenotypeName"
        placeholder="enter phenotype name"
        (ngModelChange)="addPhenotype($event)"
        [source]="phenotypeAutocompleteAdress"
        value-property-name=null
        display-property-name="name"
        path-to-data="results"
        min-chars="0"
        style="width: 100%"
        class="pt-1"/>
  </div>
    <div class="col-sm-1">
    <img src="app/images/head.png" (click)="dialog('HP:0000152')">
    </div>
    <div class="col-sm-1">
    <img src="app/images/eye.png" (click)="dialog('HP:0000478')">
    </div>
    <div class="col-sm-1">
    <img src="app/images/connect.png" (click)="dialog('HP:0003549')">
    </div>
    <div class="col-sm-1">
    <img src="app/images/bones.png" (click)="dialog('HP:0000924')">
    </div>
    <div class="col-sm-1">
    <img src="app/images/brain.png" (click)="dialog('HP:0000707')">
    </div>
    <div class="col-sm-1">
    <img src="app/images/digest.png" (click)="dialog('HP:0025031')">
    </div>
    <div class="col-sm-1">
    <img src="app/images/heart.png" (click)="dialog('HP:0001626')">
    </div>
    <div class="col-sm-1">
    <img src="app/images/muscle.png" (click)="dialog('HP:0003011')">
    </div>
  </div>
  <div *ngIf="phenotypes.length > 0" class="row pb-1 border-between">
    <div class="col-sm-5 pt-2">
      <h5>Phenotypes</h5>
      <ul>
        <li *ngFor="let phenotype of phenotypes">
          <span class="badge badge-info">{{phenotype.id}}</span> {{phenotype.name}}
          <span class="badge badge-pill badge-warning" (click)="go(phenotype.id);dialog(phenotype.id)">more</span>
          <span class="badge badge-pill badge-danger" (click)="remove(phenotype.id)">remove</span>
        </li>
      </ul>
    </div>
    <div class="col-sm-5 pt-2">
      <h5>Diseases <span class="badge badge-pill badge-default">% match</span></h5>
      <ul>
        <li *ngFor="let disease of diseases">
          {{disease.name}} <span class="badge badge-pill badge-default">{{disease.score.toFixed(1)}}</span>
        </li>
      </ul>
    </div>
    <div class="col-sm-2 pt-2">
      <button class="btn btn-default" type="button" ngxClipboard [cbContent]="genesText">copy</button>

      <h5>Genes <span class="badge badge-pill badge-default">% match</span></h5>

      <ul>
        <li *ngFor="let gene of genes">
          {{gene.name}} <span class="badge badge-pill badge-default">{{gene.score.toFixed(1)}}</span>
        </li>
      </ul>
    </div>
  </div>

  <div *ngIf="showDialog" class="dialog">
      <div *ngIf="showSpinner" class="spinner">
        <div class="bounce1"></div>
        <div class="bounce2"></div>
        <div class="bounce3"></div>
      </div>
      <h5> Phenotype </h5>
      <span class="badge badge-info">{{dialogObject.phenotype.id}}</span> {{dialogObject.phenotype.name}}
      <span class="badge badge-pill badge-primary" (click)="usePhenotype(dialogObject.phenotype.id,dialogObject.phenotype.name)">add</span>
      <div  *ngIf="dialogObject.parent != null" >
      <h5 class="mt-4"> Superclass </h5>
      <span class="badge badge-info">{{dialogObject.parent.id}}</span> {{dialogObject.parent.name}}
      <span class="badge badge-pill badge-info" (click)="dialog(dialogObject.parent.id)">more</span>
      <span class="badge badge-pill badge-primary" (click)="usePhenotype(dialogObject.parent.id,dialogObject.parent.name)">add</span>
      </div>
      <h5 *ngIf="dialogObject.children.length > 0" class="mt-4"> Subclasses </h5>
      <ul>
        <li *ngFor="let phenotype of dialogObject.children">
          <span class="badge badge-info">{{phenotype.id}}</span> {{phenotype.name}}
          <span class="badge badge-pill badge-info" (click)="dialog(phenotype.id)">more</span>
          <span class="badge badge-pill badge-primary" (click)="usePhenotype(phenotype.id, phenotype.name)">add</span>
        </li>
      </ul>
    <button (click)="close()" class="btn">Close</button>
  </div>
  <div *ngIf="showDialog" class="overlay" (click)="close()"></div>`,
styles: [`
  ng2-auto-complete, input {
    display: block; border: 0px solid #ccc; width: 300px;
  }
 .border-between > [class*='col-']:before {
 background: #999999;
 bottom: 0;
 content: " ";
 left: 0;
 position: absolute;
 width: 1px;
 top: 0;
}
.border-between > [class*='col-']:first-child:before {
 display: none;
}
.overlay {
  position: fixed;
  top: 0;
  bottom: 0;
  left: 0;
  right: 0;
  background-color: rgba(0, 0, 0, 0.5);
  z-index: 999;
}

.dialog {
  z-index: 1000;
  position: fixed;
  right: 0;
  left: 0;
  top: 20px;
  margin-right: auto;
  margin-left: auto;
  min-height: 200px;
  width: 90%;
  max-width: 820px;
  max-height: calc(100vh - 50px);
  overflow-y: auto;
  background-color: #fff;
  padding: 12px;
  box-shadow: 0 7px 8px -4px rgba(0, 0, 0, 0.2), 0 13px 19px 2px rgba(0, 0, 0, 0.14), 0 5px 24px 4px rgba(0, 0, 0, 0.12);
}

@media (min-width: 768px) {
  .dialog {
    top: 40px;
  }
}

.dialog__close-btn {
  border: 0;
  background: none;
  color: #2d2d2d;
  position: absolute;
  top: 8px;
  right: 8px;
  font-size: 1.2em;
}
`],
})
export class PhenoPanelComponent  {

  genesText : string;
  phenotypes : Phenotype[];
  diseases : Disease[];
  genes : Gene[];
  phenotypeName : string;
  phenotypeAutocompleteAdress: string = this.config.apiEndpoint + "/phenotype-autocomplete?firstLetters=:keyword&resultsCount=10";

  /* dialog */
  showDialog = false;
  showSpinner = false;
  dialogPhenotype = new BehaviorSubject<string>("");
  dialogObject : any = null;
  dialogPhenotypeOriginId : string = "";
  dialog(id : string) {
    console.log("dialog: " + id);
    this.dialogPhenotype.next(id);
  }
  close() {
    this.showDialog = false;
  }
  getPhenotype(id : string) {
    this.showSpinner = true;
    //console.log("get-phenotype-1: " + id);
    var url = this.config.apiEndpoint + "/get-phenotype-by-id?phenotypeid=" + id;
    //console.log("get-phenotype-2: " + id);
    //this.dialogObject = null;
    //console.log("get-phenotype-3: " + id);
    if (id != "") {
      //console.log("get-phenotype-4: " + id);
      this.http.get(url)
              .map(response => response.json())
              .subscribe(result =>
                { console.log("get-phenotype-subscribe: " + result.phenotype.id);
                  this.showSpinner = false;
                  this.dialogObject = result;
                  console.log(result);
                  this.showDialog = true;
                });
    }
    console.log("get-phenotype-5: " + id);
  }
  go(id : string) {
    console.log("go: " + id);
    this.dialogPhenotypeOriginId = id;
  }
  usePhenotype(id : string, name : string) {
    console.log("use-phenotype: " + id);
    this.addPhenotype(new Phenotype(id, name));
    this.remove(this.dialogPhenotypeOriginId);
    this.dialogPhenotypeOriginId = "";
    this.showSpinner = false;
    this.showDialog = false;
    this.dialogObject = null;
  }
  addPhenotype(phenotype : Phenotype) {
    console.log("add-phenotype: " + phenotype.id);
    if (phenotype.name != null) {
      this.phenotypeService.addPhenotype(phenotype);
    }
  }

  remove(id : string) {
    console.log("remove: " + id);
    for(let i = 0; i < this.phenotypes.length; i++) {
      if (this.phenotypes[i].id === id) {
        this.phenotypes.splice(i, 1);
        this.phenotypeService.phenotypes.next(this.phenotypes);
        break;
      }
    }
  }

  getGeneNames(genes : Gene[]) {
    var output : string = "";
    for (let gene of genes) {
      //console.log("name:")
      //console.log(gene.name)
      output = output + gene.name + "\n"; // 1, "string", false
    }
    this.genesText = output;
  }

  constructor(private phenotypeService: PhenotypeService,
    @Inject(APP_CONFIG) private config: IAppConfig,
    private http: Http,
    private ref: ChangeDetectorRef) {
    phenotypeService.phenotypes.subscribe(value => {
      console.log("constructorSubscribePhenotypes");
      this.phenotypes = value
    });
    phenotypeService.diseases.subscribe(value => this.diseases = value);
    phenotypeService.genes.subscribe(value => {
      this.genes = value;
      this.getGeneNames(value);
    });
    this.dialogPhenotype.subscribe(value => {
      console.log("constructorSubscribeDialogPhenotype: " + value);
      this.getPhenotype(value);
    });
  }
}
