const userModule = {
	namespaced: true,
	state: {
		userFirstName: '',
		userLastName: '',
		userName: '',
		userId: '',
		userEmail: ''
  },
	mutations: {
		setUserData (state, data) {
			state.userFirstName = data.FirstName
			state.userLastName = data.LastName
			state.userName = data.Username
			state.userId = data.UserID
			state.userEmail = data.Email
			logger.debug("Set user data for " + state.userEmail)
	  }
	},
  actions: {
    getUser() {
      logger.debug("vue.vuex.user.getUser")
			request({
				waitingText: "Verifying user credentials",
				endpoint: "user",
				method: "post",
				callback: function(data) {store.commit('user/setUserData', data)}
			})
    },
		getSamples({commit}) {
			logger.debug("vue.vuex.sample.getSamples")
			request({
				method: "get",
				waitingText: "Fetching samples information",
				endpoint: "user/privileges",
				callback: function(data){console.log(data)}
			})
		},
  } // actions:{}
}
