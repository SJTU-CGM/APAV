#!/usr/bin/perl
use strict;
use warnings;

package APAVreport;

sub getPavReport{

	my ($fa_flag, $gff_flag, $th_samples, $sample_data, $dtable, $bam_tracks) = @_;

	my $gff_track = '';
	my $gff_session = '';
	my $genome;
	if($gff_flag){
		$gff_track = 'tracks.push(getGffTrack("./browser/reference.gff.gz", "./browser/reference.gff.gz.tbi"));';
		$gff_session = 'defaultSession.view.tracks.push(getGffSession());';
	}
	if($fa_flag){
		$genome = 'const { createViewState, JBrowseLinearGenomeView } = JBrowseReactLinearGenomeView;
                        const { createElement } = React;
                        const { render } = ReactDOM;

                        let assembly = getAssembly("./browser/reference.fa.gz", "./browser/reference.fa.gz.fai", "./browser/reference.fa.gz.gzi");

                        let tracks = [];
                        '.$gff_track.'
                        '.$bam_tracks.'
                        tracks.push(getTargetTrack("./browser/target.bb"));

                        let defaultSession = getDefaultSession();
                        '.$gff_session.'

                        const state = new createViewState({
                          assembly,
                                tracks,
                                defaultSession,
                                location: initdata[3]+":"+initdata[4]+".."+initdata[5]
                        })

                        render(
                          createElement(JBrowseLinearGenomeView, { viewState: state }),
                          document.getElementById("jbrowse_linear_genome_view"),
                        )

                        pavTable.on("select", function(e, dt, type, indexes){
                                        let dataList = getPavTableSelectData(dtable, indexes);
                                        covChart.setOption(covChartOption(samples, dataList[0]));
                                        groupChart.setOption(groupChartOption(dataList[1], dataList[2]));
                                        state.session.view.navToLocString(dataList[3]+":"+dataList[4]+".."+dataList[5]);
                                })';
	}else{
		$genome = 'pavTable.on("select", function(e, dt, type, indexes){
                                        let dataList = getPavTableSelectData(dtable, indexes);
                                        covChart.setOption(covChartOption(samples, dataList[0]));
                                        groupChart.setOption(groupChartOption(dataList[1], dataList[2]));
                                })';
	}

	my $report = '<!doctype html>
<html>
<head>
<meta charset="utf-8">
<script type="text/javascript" charset="utf8" src="js/jquery-3.7.1.js"></script> 
<script type="text/javascript" charset="utf8" src="js/dataTables.js"></script> 
<script type="text/javascript" charset="utf8" src="js/dataTables.searchBuilder.js"></script> 
<script type="text/javascript" charset="utf8" src="js/searchBuilder.dataTables.js"></script> 
<script type="text/javascript" charset="utf8" src="js/dataTables.buttons.js"></script> 
<script type="text/javascript" charset="utf8" src="js/buttons.dataTables.js"></script> 
<script type="text/javascript" charset="utf8" src="js/dataTables.dateTime.min.js"></script>
<link rel="stylesheet" type="text/css" href="css/dataTables.dataTables.css">
<link rel="stylesheet" type="text/css" href="css/searchBuilder.dataTables.css">
<link rel="stylesheet" type="text/css" href="css/buttons.dataTables.css">
<link rel="stylesheet" type="text/css" href="css/dataTables.dateTime.min.css">
<script type="text/javascript" charset="utf8" src="js/dataTables.select.js"></script> 
<script type="text/javascript" charset="utf8" src="js/select.dataTables.js"></script>
<link rel="stylesheet" type="text/css" href="css/dataTables.dataTables.css">
<link rel="stylesheet" type="text/css" href="css/select.dataTables.css">
<script type="text/javascript" charset="utf8" src="js/jszip.min.js"></script> 
<script type="text/javascript" charset="utf8" src="js/pdfmake.min.js"></script> 
<script type="text/javascript" charset="utf8" src="js/vfs_fonts.js"></script> 
<script type="text/javascript" charset="utf8" src="js/buttons.html5.min.js"></script> 
<script type="text/javascript" charset="utf8" src="js/buttons.print.min.js"></script> 
<script type="text/javascript" charset="utf8" src="js/bootstrap.bundle.min.js"></script>
<link rel="stylesheet" type="text/css" href="css/bootstrap.min.css">
<script type="text/javascript" charset="utf8" src="js/echarts.js"></script> 
<script type="text/javascript" charset="utf8" src="js/report.js"></script> 
<script type="text/javascript" charset="utf8" src="js/react.production.min.js"></script> 
<script type="text/javascript" charset="utf8" src="js/react-dom.production.min.js"></script> 
<script type="text/javascript" charset="utf8" src="js/react-linear-genome-view.umd.production.min.js"></script>
<title>PAV table</title>
<style type="text/CSS">
body {
        background-color: #f4f6fa;
}
.dataTables_wrapper {
        font-size: 16px;
}
.dataTables_info {
        font-size: 14px;
        font-weight: bold
}
</style>
</head>

<body>
<br>
<div class="container">
  <div class="card">
    <div class="card-header"><b>PAV table</b></div>
    <div class="card-body">
      <div class="checkbox-group"> <b>Column selection:</b><br>
        <label><input type="checkbox" name="Chr" checked>Chr</label>
        <label><input type="checkbox" name="Start">Start</label>
        <label><input type="checkbox" name="End">End</label>
        <label><input type="checkbox" name="Length">Length</label>
        <label><input type="checkbox" name="Missing frequency" checked>Missing frequency(%)</label>
      </div>
      <table id="pav-table"  class="compact stripe">
        <thead>
          <tr>
            <th>Annnotation</th>
            <th>Chr</th>
            <th>Start</th>
            <th>End</th>
            <th>Length</th>
            <th>Missing frequency(%)</th>
            '.$th_samples.'
          </tr>
        </thead>
      </table>
    </div>
  </div>
  <br>
  <div class="row">
    <div class="col-4">
      <div class="card">
        <div class="card-header"><b>Group</b></div>
        <div class="card-body">
          <div id="group-plot" style="height:200px;"></div>
        </div>
      </div>
    </div>
    <div class="col-8">
      <div class="card">
        <div class="card-header"><b>Coverage</b></div>
        <div class="card-body">
          <div id="cov-plot" style="height:200px;"></div>
        </div>
      </div>
    </div>
  </div>
  <br>
  <div id="jbrowse_linear_genome_view"></div>
</div>
<script>

                $(document).ready(function(){

                        let samples = '.$sample_data.';
                        let dtable = '.$dtable.';

                        let covChart = echarts.init(document.getElementById("cov-plot"));
                        let groupChart = echarts.init(document.getElementById("group-plot"));
                        window.onresize = function(){
                                covChart.resize();
                                groupChart.resize();
                        }

                        let pavTable = $("#pav-table").DataTable(pavTableOption(dtable, seq(5, dtable[0].length - 1)));

			let initdata = getPavTableSelectData(dtable, [0]);
                        covChart.setOption(covChartOption(samples, initdata[0]));
                        groupChart.setOption(groupChartOption(initdata[1], initdata[2]));
                        pavTable.row(":eq(0)").select();

			'.$genome.'	
                })

</script>
</body>
</html>';
	return($report)

}

sub getSampleReport {
	my ($br_phen, $th_region, $div_phen, $idata, $sdata, $phens, $phen_init, $phen_draw) = @_;

	my $report = '<!doctype html>
<html>
<head>
<meta charset="utf-8">
<script type="text/javascript" charset="utf8" src="js/jquery-3.7.1.js"></script> 
<script type="text/javascript" charset="utf8" src="js/dataTables.js"></script> 
<script type="text/javascript" charset="utf8" src="js/dataTables.searchBuilder.js"></script> 
<script type="text/javascript" charset="utf8" src="js/searchBuilder.dataTables.js"></script> 
<script type="text/javascript" charset="utf8" src="js/dataTables.buttons.js"></script> 
<script type="text/javascript" charset="utf8" src="js/buttons.dataTables.js"></script> 
<script type="text/javascript" charset="utf8" src="js/dataTables.dateTime.min.js"></script>
<link rel="stylesheet" type="text/css" href="css/dataTables.dataTables.css">
<link rel="stylesheet" type="text/css" href="css/searchBuilder.dataTables.css">
<link rel="stylesheet" type="text/css" href="css/buttons.dataTables.css">
<link rel="stylesheet" type="text/css" href="css/dataTables.dateTime.min.css">
<script type="text/javascript" charset="utf8" src="js/dataTables.select.js"></script> 
<script type="text/javascript" charset="utf8" src="js/select.dataTables.js"></script>
<link rel="stylesheet" type="text/css" href="css/dataTables.dataTables.css">
<link rel="stylesheet" type="text/css" href="css/select.dataTables.css">
<script type="text/javascript" charset="utf8" src="js/jszip.min.js"></script> 
<script type="text/javascript" charset="utf8" src="js/pdfmake.min.js"></script> 
<script type="text/javascript" charset="utf8" src="js/vfs_fonts.js"></script> 
<script type="text/javascript" charset="utf8" src="js/buttons.html5.min.js"></script> 
<script type="text/javascript" charset="utf8" src="js/buttons.print.min.js"></script> 
<script type="text/javascript" charset="utf8" src="js/bootstrap.bundle.min.js"></script>
<link rel="stylesheet" type="text/css" href="css/bootstrap.min.css">
<script type="text/javascript" charset="utf8" src="js/echarts.js"></script> 
<script type="text/javascript" charset="utf8" src="js/report.js"></script> 
<title>Sample Information</title>
<style type="text/CSS">
body {
        background-color: #f4f6fa;
}
.dataTables_wrapper {
        font-size: 16px;
}
.dataTables_info {
        font-size: 14px;
        font-weight: bold
}
</style>
</head>

<body>
<br>
<div class="container">
  <div class="card ">
    <div class="card-header"><b>Sample table</b></div>
    <div class="card-body">
      <div class="row"> <b>Column selection:</b>
        <div class="col-4 d-flex justify-content-center align-items-center">
          <div> <b>Phenotype:</b><br>
          '.$br_phen.'
          </div>
        </div>
        <div class="col-8">
          <table id="item-table" class="compact row-border hover">
            <thead>
              <tr>
                <th></th>
                <th>Annotation</th>
                <th>Chr</th>
                <th>Start</th>
                <th>End</th>
              </tr>
            </thead>
            <tfoot>
            </tfoot>
          </table>
        </div>
      </div>
      <hr>
      <table id="sample-table" class="compact stripe" >
        <thead>
          <tr>
            <th>Sample</th>
            <th>Gender</th>
            <th>Age</th>
            <th>Location</th>
            '.$th_region.'
          </tr>
        </thead>
      </table>
    </div>
  </div>
  <br>
  <div class="row">
    <div class="col-6">
      <div class="card">
        <div class="card-header"><b>Phenotype statistics</b></div>
        <div class="card-body">
          <div id="phen">
            '.$div_phen.'
          </div>
        </div>
      </div>
    </div>
    <div class="col-6 panel">
      <div class="card">
      <div class="card-header"><b>PAV analysis result</b></div>
        <div class="card-body">
          <table id="sta-table" class="compact row-border hover">
            <thead>
              <tr>
                <th>Annotation</th>
                <th>Class</th>
                <th>Presence Number</th>
              </tr>
            </thead>
          </table>
        </div>
      </div>
    </div>
  </div>
</div>
        <script>

        $(document).ready(function(){

                let idata = '.$idata.';

                let sdata = '.$sdata.';
                let phens = '.$phens.';

                let itemTable = $("#item-table").DataTable(itemTableOption(idata));
                let sampleTable = $("#sample-table").DataTable(sampleTableOption(sdata, phens, itemTable, seq(phens.length+1, sdata[0].length-1)));

                '.$phen_init.'

                let staTable = getStaTable(sampleTable, seq(phens.length+1, sdata[0].length-1), "#sta-table", idata);

                sampleTable.on("draw", function() {
                        '.$phen_draw.'
                        staTable.destroy();
                        staTable = getStaTable(sampleTable, seq(phens.length+1, sdata[0].length-1), "#sta-table", idata);
                })

                sampleTable.draw();

        })

        </script>
</body>
</html>';

	return($report);

}

1;
