function request(params) {
  logger.debug("============== request to:")

  /* setting params */
  var waitingText = (params.waitingText === undefined) ? "Connecting with server..." : params.waitingText
  var method = (params.method === undefined) ? "get" : params.method
  var endpoint = (params.endpoint === undefined) ? "hello" : params.endpoint
  var callback = (params.callback === undefined) ? function(data){} : params.callback
  var errorCallback = (params.errorCallback === undefined) ? function(error) {
    console.log(error)
  } : params.errorCallback
  var service = (params.service === undefined) ? "API_DX" : params.service
  var reqData = (params.reqData === undefined) ? {} : params.reqData
  var headers = (params.headers === undefined) ? {
      'Authorization': 'Bearer ' + store.state.security.accessToken
  } : params.headers
  var port = 8086
  if (service == "API_DX") {port = 8086}
  if (service == "API_SECURITY") {port = 8088}

  logger.debug("" + service + "/" + endpoint)

  /* setting waiting */
  store.commit('setWaitingText', waitingText)
  store.commit('setWaitComponentVisibility', true)

  /* making request */
  axios({
    method: method,
    url: 'http://genetraps.intelliseq.pl:' + port + '/' + endpoint,
    withCredentials: true,
    crossdomain: true,
    data: Object.keys(reqData).map(function(key) {
      return encodeURIComponent(key) + '=' + encodeURIComponent(reqData[key])
    }).join('&'),
    headers: headers
  })
  .then(response => {
    store.commit('setWaitComponentVisibility', false)
    console.log(response)
    callback(response.data)
  })
  .catch(error => {
    store.commit('setWaitComponentVisibility', false)
    alert(JSON.stringify(error.response.data));
    errorCallback(error)
  })
}
