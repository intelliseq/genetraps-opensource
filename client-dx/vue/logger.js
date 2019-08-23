function Logger() {}
Logger.prototype.log = function (level, text) {
  let debug = true
  let info = true
  let warn = true
  let error = true
  if ((level == "INFO" && info) || (level == "DEBUG" && debug) || (level == "WARN" && warn) || (level == "ERROR" && error)){
    today = new Date();
    dd = today.getDate();
    mm = today.getMonth()+1; //January is 0!
    yyyy = today.getFullYear();
    if (dd < 10) {dd = '0' + dd;}
    if (mm < 10) {mm = '0' + mm;}
    var date_string = dd + '/' + mm + '/' + yyyy + " " + today.getHours() + ":" + today.getMinutes()
    console.log("LOG: " + date_string + " " + level + " " + text)
  }
}
Logger.prototype.debug = function(text) {
  this.log("DEBUG", text)
}
Logger.prototype.info = function(text) {
  this.log("INFO", text)
}
Logger.prototype.warn = function(text) {
  this.log("WARN", text)
}
Logger.prototype.error = function(text) {
  this.log("ERROR", text)
}
/*Logger.prototype.log = function(text) {
  this.debug(text)
}*/
logger = new Logger()
