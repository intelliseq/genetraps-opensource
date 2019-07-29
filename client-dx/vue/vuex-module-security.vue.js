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
	        "username": credentials.login,
	        "password": credentials.password
	      },
				headers: {
          'authorization': 'Basic d2ViX2FwcDpzZWNyZXQ=',
          'content-Type': 'application/x-www-form-urlencoded'
      	},
				callback: function(data) {
					logger.debug("vue.app.loginWithRefreshToken succesfull login")
					store.dispatch('security/login', data)
					router.push("/")
				} // calback:
			}) // request()
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
					store.dispatch("security/login", data)
				}, // calback:
				errorCallback: function(error) {
					logger.debug("vue.app.loginWithRefreshToken succesfull failure")
					router.push("/login")
				}
			}) // request()
    } // loginWithRefreshToken()

  } // actions:
} // securityModule:
