const samplesModule = {
	namespaced: true,
	state: {
		samples: {}
  },
	mutations: {
		setSamples (state, samplesPrivileges) {
			logger.debug("vue.vuex.samples.setSamples")
			var samples = []
			for (var key in samplesPrivileges) {
    		if (samplesPrivileges.hasOwnProperty(key)) {
					samples.push({id: key, access: samplesPrivileges[key]})
    		}
			}
			state.samples = samples
	  }
	},
  actions: {
		getSamples({commit}) {
			logger.debug("vue.vuex.sample.getSamples")
			request({
				waitingText: "Fetching samples information",
				endpoint: "user/privileges",
				callback: function(data){
					store.commit('samples/setSamples', data)
				}
			})
		},
		showSamples({commit}) {
			logger.debug("vue.vuex.sample.showSamples")
			store.dispatch("samples/getSamples")
			store.commit("setToolbarVisibility", true)
			router.push("/samples")
		},
	}
}
