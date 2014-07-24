package IT_config_utility;

use Exporter;        # load Exporter module

@ISA=qw(Exporter);   # Inherit from Exporter
@EXPORT=qw(appendAttributes loadConfig);

sub appendAttributes($$)
{
        my ($file, $config) = @_;

        foreach my $key (sort {$a cmp $b} keys %{$config->{"ATTRIBUTE"}})
        {
                my $value = $config->{"ATTRIBUTE"}->{$key};
                $file->print("data\$".$key." <- \"$value\"\n");
        }
}

sub loadConfig($$)
{
        my ($file_name, $config) = @_;

        my $file = IO::File->new($file_name, "r");

        while(my $line = <$file>)
        {
                chomp($line);

                my ($name, $value) = split(/\t/, $line);
                my ($att, $att_val) = split(":", $value);

                $config->{$name}->{$att} = $att_val;
        }
}



1;
