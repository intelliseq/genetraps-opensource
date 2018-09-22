const securityModule = {
	namespaced: true,
	state: {
      client_id: "web_app",
			accessToken: '',
    },
		mutations: {
				updateToken: (state, token) => {
					logger("DEBUG", "token updated")
					state.accessToken = token
				},
		},
    actions: {
			login({commit}, response) {
				logger("DEBUG", "vue.vuex.security.login")
				Vue.cookies.set("refresh_token", response.refresh_token, "7D")
				store.commit('security/updateToken', response.access_token)
				store.dispatch('user/getUser')
				//app.$router.push("/")
			},
      loginWithCredentials({commit}, credentials) {
				store.commit('setWaitingText', "Signing in")
				store.commit('setWaitingVisibility', true)
        logger("DEBUG", "vue.vuex.security.loginWithCredentials")
        var reqData = {
          "grant_type": "password",
          "client_id": store.state.security.client_id,
          "username": credentials.login,
          "password": credentials.password
        }
        axios({
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
          logger("DEBUG",  "vue.login.getToken")
          logger("DEBUG",  "succesfull login")
					store.commit('setWaitingVisibility', false)
          store.dispatch('security/login', response.data)
          //console.log(response.data.access_token)
        })
        .catch(e => {
          console.log(e)
					store.commit('setWaitingVisibility', false)
          //this.error = true
          //this.dialog_visibility = false
          //this.error_message = "Error: " + e.response.data.error
          //this.error_description = "Error description: " + e.response.data.error_description
          //this.alert_visibility = true
        })
      },
      loginWithRefreshToken({commit}) {
				store.commit('setWaitingText', "Signing in")
				store.commit('setWaitingVisibility', true)
        //this.$hub.$emit('start_waiting', "Signing in")
        var reqData = {
          "grant_type": "refresh_token",
          "client_id": "web_app",
          "refresh_token": Vue.cookies.get("refresh_token")
        }
				//console.log(reqData)
        axios({
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
          logger("DEBUG", "vue.app.loginWithRefreshToken succesfull login")
					store.commit('setWaitingVisibility', false)
					store.dispatch('security/login', response.data)
          //this.$store.commit('login', response.data);
          //console.log(response.data.access_token)
        })
        .catch(e => {
          logger("ERROR", e)
					store.commit('setWaitingVisibility', false)
          //this.error = true
          //this.dialog_visibility = false
          //this.error_message = "Error: " + e.response.data.error
          //this.error_description = "Error description: " + e.response.data.error_description
          //this.alert_visibility = true
        })
      }

    }
  }
