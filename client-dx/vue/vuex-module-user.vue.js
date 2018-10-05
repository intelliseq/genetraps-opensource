const userModule = {
	namespaced: true,
	state: {
        userFirstName: '',
        userLastName: '',
        userName: '',
        userId: '',
        userEmail: ''
    },
    actions: {
      getUser() {
        logger("DEBUG", "vue.vuex.user.getUser")
				logger("DEBUG", "TOKEN: " + store.state.security.accessToken)
				store.commit('setWaitingText', "Verifying user credentials")
				store.commit('setWaitingVisibility', true)
	      axios({
          method: 'post', //you can set what request you want to be
	        url: 'http://genetraps.intelliseq.pl:8086/user',
	        withCredentials: true,
	        crossdomain: true,
	        headers: {
	            'Authorization': 'Bearer ' + store.state.security.accessToken},
	        })
	        .then(response => {
	          console.log(response)
	          //logger("DEBUG",  "vue.login.getToken")
	          //logger("DEBUG",  "succesfull login")
						store.commit('setWaitingVisibility', false)
	          //store.dispatch('security/login', response.data)
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
      }

  }
