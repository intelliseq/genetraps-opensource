const leftToolbarComponent = {
  template:
`
<v-flex v-if="toolbarVisibility" d-flex xs12 sm2 md2>
<v-card color="teal lighten-3" >
  <v-layout column px-2 py-0>

    <toolbar-button-component
    v-for="(item, index) in todos"
    v-bind:item="item"
    v-bind:index="index"
    v-bind:key="item.id"
  ></toolbar-button-component>
  </v-layout>
</v-card>
</v-flex>
`,
components: {
'toolbar-button-component' :toolbarButtonComponent,
},
    computed: {
      ...Vuex.mapState(['toolbarVisibility'])
    },
    //  Vuex.mapState('user', ['userEmail'])
    //}
    data: function() {return {

    todos: [
      {
        id: 1,
        title: 'CREATE',
        icon: 'add_circle',
        tooltip: 'create new sample',
      },
      {
        id: 2,
        title: 'Two',
        icon: 'add_circle',
        tooltip: 'tooltip',
      }
    ],

    }
  },
    created: function () {
      logger.debug("vue.toolbar-left.created")
    }
}
