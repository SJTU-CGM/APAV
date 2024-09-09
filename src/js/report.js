
function sum(arr){
	return arr.reduce(function(prev, curr, idx, arr){
		return Number(prev) + Number(curr);
	})
}

function seq(start, end){
	let arr = [];
	for (let i = start; i <= end; i++) {
	  arr.push(i);
	}
	return arr;
}

function colCheckbox(api, name){
			$("input[name='"+name+"']").on('keyup click', function(){
					if($("input[name='"+name+"']").prop("checked")){
					   api.column( ":contains("+name+")" ).visible( true );	
					} else {
						api.column( ":contains("+name+")").visible( false );				 
					}
				})
		}

function getAssembly(fa, fai, gzi){
	return {
			  name: 'Pangenome',
			  sequence: {
				type: 'ReferenceSequenceTrack',
				trackId: 'Pan-ReferenceSequenceTrack',
				adapter: {
				  type: 'BgzipFastaAdapter',
				  fastaLocation: {
					uri: fa,
					locationType: 'UriLocation',
				  },
				  faiLocation: {
					uri: fai,
					locationType: 'UriLocation',
				  },
				  gziLocation: {
					uri: gzi,
					locationType: 'UriLocation',
				  },
				},
			  },
			}
}

function getGffTrack(gff, tbi){
	return {
    type: 'FeatureTrack',
    trackId:
      'Pan-AnnotationTrack',
    name: 'Annotation',
    assemblyNames: ['Pangenome'],
    adapter: {
      type: 'Gff3TabixAdapter',
      gffGzLocation: {
        uri: gff,
        locationType: 'UriLocation',
      },
      index: {
        location: {
          uri: tbi,
          locationType: 'UriLocation',
        },
        indexType: 'TBI',
      },
    },
    renderer: {
      type: 'SvgFeatureRenderer',
    },
  }
}


function getBamTrack(name, bam, bai){
	return {
  trackId: "pan-AlignmentTrack-"+name,
  name: "Alignment-"+name,
  assemblyNames: ['Pangenome'],
  type: "AlignmentsTrack",
  adapter: {
    "type": "BamAdapter",
    "bamLocation": {
      "uri": bam
    },
    index: {
      "location": {
        "uri": bai
      }
    }
  }
}
}

function getTargetTrack(bb){
	return {
    "type": "FeatureTrack",
    "trackId": "Pan-TargetTrack",
    "name": "Target region",
    "adapter": {
      "type": "BigBedAdapter",
      "bigBedLocation": {
        "uri": bb,
        "locationType": "UriLocation"
      }
    },
    "assemblyNames": [
      "Pangenome"
    ]
  }
}

function getDefaultSession(){
	return {
          name: 'my session',
          view: {
            type: 'LinearGenomeView',
            tracks: [
              {
                type: 'ReferenceSequenceTrack',
                configuration: 'Pan-ReferenceSequenceTrack',
                displays: [
                  {
                    type: 'LinearReferenceSequenceDisplay',
                    configuration:
                      'Pan-ReferenceSequenceTrack-LinearReferenceSequenceDisplay',
                  },
                ],
              },
				{
				type: 'FeatureTrack',
				configuration: 'Pan-TargetTrack',
				displays: [
				  {
					type: 'LinearBasicDisplay',
					configuration: 'Pan-TargetTrack-LinearBasicDisplay',
					  height:30
				  },
				],
			  },
            ],
          },
        }
}

function getGffSession(){
	return {
				type: 'FeatureTrack',
				configuration: 'Pan-AnnotationTrack',
				displays: [
				  {
					type: 'LinearBasicDisplay',
					configuration: 'Pan-AnnotationTrack-LinearBasicDisplay',
				  },
				],
			  }
}

function pavTableOption(dtable, idxs){
	return {
					data:dtable,
					dom:'Qrtip',
					paging: false,
                    scrollX: true,
                    scrollY: 200,
					columnDefs:[
						{
							targets: idxs,
							render: function(data, type, row){
								let arr = data.split("/");
								return arr[0];
							}
						},
						{
							targets: [2,3,4],
							visible: false
						}
					],
					language:{
						searchBuilder:{
							title:{
								0:'<b>Advanced Filter:</b>',
								_:'<b>Advanced Filter:</b>'
							}
						}
					},
					select: true,
					order: [],
					initComplete: function () {		
						let cols = ["Chr", "Start", "End", "Length", "Missing frequency"]
						var api = this.api();	
						cols.forEach(function(item){
							colCheckbox(api, item)
						})
					}
				}
}

function getPavTableSelectData(dtable, indexes){
	let cur = dtable[indexes[0]];
	let curd = cur.toSpliced(0,6);
	let starts = cur[2].split(",");
	let ends = cur[3].split(",");
	
					let sdata = [];
					let pdata = [];
					let adata = [];
					curd.forEach(function(i){
						let tmp = i.split("/");
						let c = tmp[0] == 1 ? 'steelblue' : '#F08080';	
						tmp[0] == 1 ? pdata.push(tmp[1]) : adata.push(tmp[1]);
						sdata.push({value:tmp[1],
								   itemStyle:{
									   color: c
								   }})
					})
	
	
	return [sdata, adata, pdata, cur[1], starts[0], ends[ends.length-1]]
}

function covChartOption(samples,sdata){
	return {
						xAxis:{
							type: 'category',
							data: samples,
							axisTick:{
								alignWithLabel:true,
								interval: 0
							},
							axisLabel:{
								rotate:45
							}
						},
						yAxis:{
							type: 'value'
						},
						series: [
							{
								data: sdata,
								type: 'bar'
							}
						]
					}
}


function groupChartOption(adata, pdata){
	let groupX = ["Absence", "Presence"];
	return {
							  dataset: [
								{
								  // prettier-ignore
								  source: [
											adata,pdata
										]
								},
								{
								  transform: {
									type: 'boxplot',
									  config:{itemNameFormatter:function(params){
										  return groupX[params.value]
									  }
											 }
								  }
								},
								{
								  fromDatasetIndex: 1,
								  fromTransformResult: 1
								}
							  ],
							  tooltip: {
								trigger: 'item',
								axisPointer: {
								  type: 'shadow'
								}
							  },
							  grid: {
								left: '10%',
								right: '10%',
								bottom: '15%'
							  },
							  xAxis: {
								type: 'category',
								  axisTick:{
								alignWithLabel:true,
								interval: 0
							}
							  },
							  yAxis: {
								type: 'value',
								name: 'Coverage'
							  },
							  series: [
								{
								  name: 'boxplot',
								  type: 'boxplot',
									colorBy:'data',
									color:["#F08080", "steelblue"],
								  datasetIndex: 1
								},
								{
								  name: 'outlier',
								  type: 'scatter',
									itemStyle:{color:'gray'},
									symbolSize:5,
								  datasetIndex: 2
								}
							  ]
					
					}
}



function itemTableOption(d){
	return {
			data: d,
			dom:'ft',
			paging: false,
			scrollY: 100,
			columnDefs: [ {
				orderable: false,
				className: 'select-checkbox',
				targets:   0
			} ],
			select: {
				style:    'multi',
				selector: 'td:first-child'
			}
		}
}

function sampleTableOption(sdata, phens, itemTable, idxs){
	return {
					data:sdata,
					dom:'Qrtip',
					paging: false,
                    scrollY: 100,
					columnDefs:[
						{
							targets: idxs,
							visible : false
						}
					],
					language:{
						searchBuilder:{
							title:{
								0:'<b>Advanced Filter:</b>',
								_:'<b>Advanced Filter:</b>'
							}
						}
					},
					select: true,
					initComplete: function () {
						let api = this.api();	
						phens.forEach(function(item){
							colCheckbox(api,item);
						})
						
						itemTable
							.on( 'select', function ( e, dt, type, indexes ) {
							let rowData = itemTable.rows( indexes ).data().toArray();
							api.column( ':contains("'+rowData[0][1]+'")' ).visible( true );	
						} )
						.on( 'deselect', function ( e, dt, type, indexes ) {
							let rowData = itemTable.rows( indexes ).data().toArray();
							api.column( ':contains("'+rowData[0][1]+'")' ).visible( false );	
						} )
					}
				}
}


function getPhenData(table, colidx) {
    		let counts = {}
 
    		table
        		.column(colidx, { search: 'applied' })
        		.data()
        		.each(function (val) {
            		if (counts[val]) {
                		counts[val] += 1;
					} else {
                		counts[val] = 1;
            		}
        		});
 
			let counts2 = {}		
			Object.keys(counts).sort().forEach(function(key){
				counts2[key] = counts[key];
			});
			
    		return $.map(counts2, function (val, key) {
        		return {
            		name: key,
            		value: val,
        		};
    		});
		}

function phenChartOption(d) {
	return {
			 tooltip: {
				trigger: 'item',
							 show: false
			  },
			  legend: {
				top: '5%',
				left: 'center'
			  },
			  series: [
				{
				  type: 'pie',
				  radius: ['40%', '70%'],
				  itemStyle: {
					borderRadius: 10,
					borderColor: '#fff',
					borderWidth: 2
				  },
				  label: {
					show: false,
					position: 'center'
				  },
				  emphasis: {
					label: {
					  show: true,
					  fontSize: 20,
					  fontWeight: 'bold',
						formatter(param) {
										return param.name + '\n(' + param.percent  + '%)';
									}
					}
				  },
				  data: d
				}
			  ]
		}
	
}


function getStaTable(source, idxs, table, info){
			let sta = [];
			source.columns(idxs, { search: 'applied' }).data().each(function (val, idx) {
				if(val.length > 0){
					sta.push([info[idx][1], sum(val) == val.length ? "Core" : "Distributed", sum(val)])
				}				
			})
			
			return $(table).DataTable({
				data: sta,
				layout: {
        topStart: {
            buttons: ['searchBuilder']
        },
					topEnd: null,
					bottomStart: null,
					bottomEnd:{
						buttons: ['excel', 'csv']
					}
    },
				paging: false,
				scrollY: 100,
				order:[[2,'asc']]
			})
			
		}
