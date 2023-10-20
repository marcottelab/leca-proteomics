use strict;
use warnings;
use LWP::UserAgent;
use HTTP::Date;

my $slc_prefix = $ARGV[0];  # include descriptive names for better file naming
my $slc_id = $ARGV[1]; # uniprot identifier for subcellular localization

my $query = "https://rest.uniprot.org/uniprot/stream?query=$slc_id+reviewed:true&fields=accession,id,organism_id,protein_families,go,go_p,go_c,go_f,go_id,cc_subcellular_location,lineage,lineage_ids&format=tsv";

my $file = $slc_prefix . '_' . $slc_id . '.tsv';

my $contact = 'rachaelcox@utexas.edu'; # Please set a contact email address here to help us debug in case of problems (see https://www.uniprot.org/help/privacy).
my $agent = LWP::UserAgent->new( agent => "libwww-perl $contact" );
my $response = $agent->mirror( $query, $file );

if ( $response->is_success ) {
	my $results      = $response->header('X-Total-Results');
	my $release      = $response->header('X-UniProt-Release');
	my $release_date = $response->header('X-UniProt-Release-Date');
	print
"Downloaded metadata for SLC ID: $slc_id from UniProt release $release ($release_date) to file $file\n";
}
elsif ( $response->code == HTTP::Status::RC_NOT_MODIFIED ) {
print "Data for subcellular term $slc_id is up-to-date.\n";
}
else {
	die 'Failed, got '
	. $response->status_line . ' for '
	. $response->request->uri . "\n";
}
