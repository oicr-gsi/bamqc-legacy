requires 'JSON', '2.00'; # 2.0 or newer

on 'test' => sub {
    requires 'Text::Diff';
    requires 'File::Slurp';
};

requires 'SeqWare::Html' => 0.01, git => 'git@github.com:oicr-gsi/gsi-website.git', rev => '53f543e';
