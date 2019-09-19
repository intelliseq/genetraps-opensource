const sampleComponent = {
  template:`
  <v-layout>
  <v-flex ma-1>
  <v-card>
    <v-list two-line subheader>
      <v-subheader inset>Sample {{$route.params.id}}</v-subheader>
        <v-list-tile v-for="(value, key) in sample.properties" :key="key" avatar @click="">
            <v-list-tile-avatar>
              <v-icon>label</v-icon>
            </v-list-tile-avatar>
            <v-list-tile-content>
              <v-list-tile-title>{{ key }}</v-list-tile-title>
              <v-list-tile-sub-title>{{ value }}</v-list-tile-sub-title>
            </v-list-tile-content>
            <v-list-tile-action>
              <v-btn icon>
                <v-icon color="grey lighten-1">info</v-icon>
              </v-btn>
            </v-list-tile-action>
          </v-list-tile>


<!--           <v-list-tile v-for="file in sample.files" avatar @click="">
                      <v-list-tile-avatar>
                        <v-icon>label</v-icon>
                      </v-list-tile-avatar>
                      <v-list-tile-content>
                        <v-list-tile-title>{{ file.name }}</v-list-tile-title>
                        <v-list-tile-sub-title>file</v-list-tile-sub-title>
                      </v-list-tile-content>
                      <v-list-tile-action>
                        <v-btn icon>
                          <v-icon color="grey lighten-1">info</v-icon>
                        </v-btn>
                      </v-list-tile-action>
                    </v-list-tile> -->


          <v-list-tile>
          <v-list-tile-avatar>
            <v-icon color="green darken-2">add_circle</v-icon>
          </v-list-tile-avatar>
          <v-text-field @keyup.enter="add" v-model="newKey" label="Property Name"></v-text-field></td>
          <v-text-field @keyup.enter="add" v-model="newValue" label="Property Value"></v-text-field>
          </v-list-tile>


           <v-list-tile>
          <v-list-tile-avatar><v-icon color="green darken-2">attach_file</v-icon> </v-list-tile-avatar>
          <input type="file" multiple @change="handleFilesUpload"></input>
          <button @click="submitFiles">Submit</button>
         </v-list-tile>




        </v-list>
        <v-list>
        <div v-for="file in files">
          PLIK:{{file.name}}
        </div>

        </v-list>
  </v-card>
  </v-flex>
  </v-layout>

    `,
    data: function () {
      return {
        newKey: undefined,
        newValue: undefined,
        files: []
      }
    },
    computed: {
      ...Vuex.mapState('sample', ['sample'])
    },
    methods: {
      add(){
    	   Vue.set(this.sample.properties, this.newKey, this.newValue)
         this.newName = undefined;
         this.newNumber = undefined;
         console.log(this.sample)
      },

       submitFiles(){

             let formData = new FormData();

               for( var i = 0; i < this.files.length; i++ ){
                 let file = this.files[i];

                 formData.append('files[' + i + ']', file);
                 //Vue.set(this.sample.files, "file", file)

                 //Vue.set(this.sample.files, i, file)
               }

               /*
                 Make the request to the POST /multiple-files URL
               */

//                axios.post( '/sample/{id}/file/upload',
//                  formData,
//                  {
//                    headers: {
//                        'Content-Type': 'multipart/form-data'
//                    }
//                  }
//                ).then(function(){
//                  console.log('SUCCESS!!');
//                })
//                .catch(function(){
//                  console.log('FAILURE!!');
//                });

             },

             handleFilesUpload(event){
               var files = []
               for (i = 0; i < event.target.files.length; i++) {
                 var file = {}
                 file["name"] = event.target.files[i].name
                 file["file"] = event.target.files[i];
                 files.push(file)
               }
               this.files = files
               //this.sample.files = files;
               //console.log(event.target.files)
               //console.log(event.target.files[0])
               console.log("success")
             }

    },
    created: function () {
      logger.debug("vue.sample.created")
      logger.debug(this.currentSampleId)
    }
}
