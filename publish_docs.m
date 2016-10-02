% Script to automate generating documentation for demo cases. Output
% locations and filenames are chosed to play well with Github Pages.

publish('even_stream_demo.m', 'format', 'html', 'outputDir', 'docs');
movefile('docs/even_stream_demo.html', 'docs/index.html');