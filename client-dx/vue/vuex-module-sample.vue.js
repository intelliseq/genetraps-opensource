const sampleModule = {
	namespaced: true,
	state: {
		sampleId: ''
  },
	mutations: {
	},
  actions: {
		createSample({commit}) {
			logger.debug("vue.vuex.sample.createSample")
			request({
				waitingText: "Creating sample",
				endpoint: "mkdir",
				callback: function(data){console.log(data)}
			})
		},
		getSamples({commit}) {
			logger.debug("vue.vuex.sample.getSamples")
			request({
				waitingText: "Fetching samples information",
				endpoint: "user/samples",
				callback: function(data){console.log(data)}
			})
		},
	}
}
