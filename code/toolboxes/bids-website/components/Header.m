function html = Header(props)

html = '';
html = [html, '\n', '<html lang="en" data-bs-theme="light">'];
html = [html, '\n', '<head>'];
html = [html, '\n', '   <meta charset="utf-8">'];
html = [html, '\n', '   <meta name="viewport" content="width=device-width, initial-scale=1">'];
html = [html, '\n', '   <meta name="author" content="Rick Wassing, Woolcock Institute of Medical Research, Sydney, Autralia">'];
html = [html, '\n', '   <title>', props.title, '</title>'];
html = [html, '\n', '   <link rel="canonical" href="https://getbootstrap.com/">'];
html = [html, '\n', '   <!-- jQuery -->'];
html = [html, '\n', '   <link href="https://cdn.datatables.net/1.13.1/css/jquery.dataTables.min.css" rel="stylesheet">'];
html = [html, '\n', '   <script   src="https://code.jquery.com/jquery-3.6.2.min.js"   integrity="sha256-2krYZKh//PcchRtd+H+VyyQoZ/e3EcrkxhM8ycwASPA="   crossorigin="anonymous"></script>'];
html = [html, '\n', '   <script type="text/javascript" src="https://cdn.datatables.net/1.13.1/js/jquery.dataTables.min.js"></script>'];
html = [html, '\n', '   <!-- Bootstrap core CSS -->'];
html = [html, '\n', '   <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-rbsA2VBKQhggwzxH7pPCaAqO46MgnOM80zW1RWuH61DGLwZJEdK2Kadq2F9CUG65" crossorigin="anonymous">'];
html = [html, '\n', '   <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/js/bootstrap.bundle.min.js" integrity="sha384-kenU1KFdBIe4zVF0s0G1M5b4hcpxyD9F7jL+jjXkk+Q2h455rYXK/7HAuoJl+0I4" crossorigin="anonymous"></script>'];
html = [html, '\n', '   <!-- Styling -->'];
html = [html, '\n', '   <style>'];
html = [html, '\n', '       body {'];
html = [html, '\n', '           font-size: 10pt;'];
html = [html, '\n', '           padding: 12px;'];
html = [html, '\n', '           height: 100%%;'];
html = [html, '\n', '           background-color: 100%%;'];
html = [html, '\n', '           background-attachment: fixed;'];
html = [html, '\n', '           background-repeat: no-repeat;'];
html = [html, '\n', '           background-image: '];
html = [html, '\n', '               linear-gradient(180deg, rgba(248, 249, 250, 0.01), rgba(248, 249, 250, 1) 85%%),'];
html = [html, '\n', '               radial-gradient(ellipse at top left, rgba(13, 110, 253, 0.5), transparent 50%%),'];
html = [html, '\n', '               radial-gradient(ellipse at top right, rgba(255, 228, 132, 0.5), transparent 50%%),'];
html = [html, '\n', '               radial-gradient(ellipse at center right, rgba(113, 44, 249, 0.5), transparent 50%%),'];
html = [html, '\n', '               radial-gradient(ellipse at center left, rgba(214, 51, 132, 0.5), transparent 50%%);'];
html = [html, '\n', '       }'];
html = [html, '\n', '       h1 { font-size: 20pt; }'];
html = [html, '\n', '       h2 { font-size: 18pt; }'];
html = [html, '\n', '       h3 { font-size: 16pt; }'];
html = [html, '\n', '       h4 { font-size: 14pt; }'];
html = [html, '\n', '       h5 { font-size: 12pt; }'];
html = [html, '\n', '       h6 { font-size: 10pt; }'];
html = [html, '\n', '       .w-8 { max-width: 8.333333%%; }'];
html = [html, '\n', '       .w-17 { max-width: 16.666666%%; }'];
html = [html, '\n', '       .w-33 { max-width: 33.333333%%; }'];
html = [html, '\n', '   </style>'];
html = [html, '\n', '</head>'];
html = [html, '\n', '<body>'];
html = [html, '\n', '   <main>'];
html = [html, '\n', '       <div class="container-fluid">'];

end