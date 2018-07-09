"use strict";
var Disease = (function () {
    //public id : string;
    /*public getId() : string {
      return this.id;
    }*/
    function Disease(name, match) {
        this.name = name;
        this.match = match;
    }
    return Disease;
}());
exports.Disease = Disease;
//# sourceMappingURL=disease.js.map