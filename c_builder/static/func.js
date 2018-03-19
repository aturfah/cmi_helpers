function select_dir() {
    $.get("/files/",
        success = function (data) {
            if (!data['status']) {
                alert("Load aborted.")
                return;
            }

            sample_data = data['sample']
            load_sample(sample_data)
        })
}

function load_C() {
    $.get("/c/", 
        success = function (data) {
            if (!data['status']) {
                alert("Load aborted.")
                return;
            }

            sample_data = data['sample']
            load_sample(sample_data)
        })
}

function load_sample(sample_data) {
    // $('#sample_dir').html('Loaded');
    var sample_table = $('#sample_table')
    sample_table.empty()

    var header_row = $('<tr>')
    header_row.append('<td>Regular Expression</td>')
    header_row.append('<td>Reg</td>')
    header_row.append('<td>Modality</td>')
    header_row.append('<td>Label</td>')
    header_row.append('<td>Incremental</td>')
    header_row.append('<td>Ignore</td>')
    
    sample_table.append(header_row)
    populate_form(sample_table, sample_data, "")

}

function populate_form(form_table, data, base_name) {
    for (file_ in data) {
        // console.log(data[file_])
        if (file_ === "name") {
            continue
        }
        
        if (data[file_].hasOwnProperty("file") && data[file_]["file"] === 1) {
            var tr = $('<tr>');
            // Regular Expression
            tr.append('<td><input name="regex" type="text" size=35 value="REPLACE"></td>'.replace("REPLACE", base_name + "/" + file_))
            
            // Reg
            var reg_val = ""
            if ("reg" in data[file_]) {
                reg_val = data[file_]["reg"]
            }
            tr.append('<td><input name="reg" type="text" value="REPLACE"></td>'.replace("REPLACE", reg_val))

            // Modality
            var mod_val = ""
            if ("mod" in data[file_]) {
                mod_val = data[file_]["mod"]
            }
            tr.append('<td><input name="modality" type="text" value="REPLACE"></td>'.replace("REPLACE", mod_val))

            // Label
            var lab_val = ""
            if ("label" in data[file_]) {
                mod_val = data[file_]["label"]
            }
            tr.append('<td><input name="label" type="text" value="REPLACE"></td>'.replace("REPLACE", lab_val))

            // Incremental
            var incr_val = ""
            if ("incremental" in data[file_]) {
                incr_val = data[file_]["incremental"]
            }
            tr.append('<td><input name="incremental" type="checkbox" value="REPLACE"></td>'.replace("REPLACE", incr_val))

            // Ignore
            tr.append('<td><input name="ignore" type="checkbox"></td>')

            form_table.append(tr)
        } else {
            populate_form(form_table, data[file_], "/" + file_)
        }
    }
}

function export_file() {
    console.log("Exporting file...")
    var post_data = {}

    var data_table = $("#sample_table")

    data_table.children().each(function(i, row) {
        if(i === 0) {
            // Header row
            return
        }
        var datum = {}
        var row = $(row)

        datum['regex'] = row.find('input[name=regex]').val();
        datum['reg'] = row.find('input[name=reg]').val();
        datum['label'] = row.find('input[name=label]').val();
        datum['modality'] = row.find('input[name=modality]').val()
        datum['incremental'] = row.find('input[name=incremental]').is(":checked");
        datum['ignore'] = row.find('input[name=ignore]').is(":checked");
        
        var regex = datum['regex']
        post_data[regex] = datum
    });

    console.log(post_data)

    $.post("/export/",
        data = JSON.stringify(post_data),
        success = function (data) {
            if (!data['status']) {
                alert("Export failed.")
                return;
            }

           console.log("Success!")
        })
}

function clear_dir() {
    var sample_table = $('#sample_table')
    sample_table.empty()

    $.get("/clear/")
}