"use strict";
var Gene = (function () {
    //public id : string;
    /*public getId() : string {
      return this.id;
    }*/
    function Gene(name, match) {
        this.name = name;
        this.match = match;
    }
    return Gene;
}());
exports.Gene = Gene;
//# sourceMappingURL=gene.js.map