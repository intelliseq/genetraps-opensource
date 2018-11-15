const securityModule = {
	namespaced: true,
	state: {
    client_id: "web_app",
		accessToken: '',
  },
	mutations: {
		updateToken: (state, token) => {
			logger.debug("token updated")
			state.accessToken = token
		},
	},
  actions: {
		login({commit}, response) {
			logger.debug("vue.vuex.security.login")
			Vue.cookies.set("refresh_token", response.refresh_token, "7D")
			store.commit('security/updateToken', response.access_token)
			store.dispatch('user/getUser')
		},
    loginWithCredentials({commit}, credentials) {
      logger.debug("vue.vuex.security.loginWithCredentials")
			request({
				waitingText: "Signing in",
				service: "API_SECURITY",
				endpoint: "oauth/token",
				method: "post",
				reqData: {
	        "grant_type": "password",
	        "client_id": store.state.security.client_id,
	        "username": credentials.login,
	        "password": credentials.password
	      },
				headers: {
          'Authorization': 'Basic d2ViX2FwcDpzZWNyZXQ=',
          'Content-Type': 'application/x-www-form-urlencoded'
      	},
				callback: function(data) {
					logger.debug("vue.app.loginWithRefreshToken succesfull login")
					store.dispatch('security/login', data)
					router.push("/")
				} // calback:
			}) // request()
      /*var reqData = {
        "grant_type": "password",
        "client_id": store.state.security.client_id,
        "username": credentials.login,
        "password": credentials.password
      }*/
      /*axios({
        method: 'post', //you can set what request you want to be
        url: 'http://genetraps.intelliseq.pl:8088/oauth/token',
        withCredentials: true,
        crossdomain: true,
        data: Object.keys(reqData).map(function(key) {
          return encodeURIComponent(key) + '=' + encodeURIComponent(reqData[key])
        }).join('&'),
        headers: {
          'Authorization': 'Basic d2ViX2FwcDpzZWNyZXQ=',
          'Content-Type': 'application/x-www-form-urlencoded' },
      })
      .then(response => {
        logger.debug( "vue.login.getToken")
        logger.debug( "succesfull login")
				store.commit('setWaitingVisibility', false)
        store.dispatch('security/login', response.data)
        //console.log(response.data.access_token)
      })
      .catch(e => {
        console.log(e)
				store.commit('setWaitingVisibility', false)
      })*/
    },
    loginWithRefreshToken({commit}) {
			request({
				waitingText: "Signing in",
				service: "API_SECURITY",
				endpoint: "oauth/token",
				method: "post",
				reqData: {
          "grant_type": "refresh_token",
          "client_id": "web_app",
          "refresh_token": Vue.cookies.get("refresh_token")
        },
				headers: {
          'Authorization': 'Basic d2ViX2FwcDpzZWNyZXQ=',
          'Content-Type': 'application/x-www-form-urlencoded'
      	},
				callback: function(data) {
					logger.debug("vue.app.loginWithRefreshToken succesfull login")
					store.dispatch('security/login', data)
				} // calback:
			}) // request()
    } // loginWithRefreshToken()
  } // actions:
} // securityModule:
