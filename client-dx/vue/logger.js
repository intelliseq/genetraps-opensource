function logger(level, text) {
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
