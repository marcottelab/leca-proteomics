use strict;
use warnings;
use LWP::UserAgent;
use HTTP::Date;

my $organism_id = $ARGV[0]; # uniprot identifier of organism
my $tax_id = $ARGV[1]; # ncbi identifier of organism

#my $query = "https://rest.uniprot.org/proteomes/stream?query=$organism_id+$tax_id&fields=upid,organism,organism_id,protein_count,mnemonic,lineage&format=tsv";
my $query = "https://rest.uniprot.org/proteomes/stream?query=$organism_id+$tax_id+reference:true&fields=upid,organism,organism_id,protein_count,mnemonic,lineage&format=tsv";

my $file = $organism_id . '.tsv';

my $contact = 'rachaelcox@utexas.edu'; # Please set a contact email address here to help us debug in case of problems (see https://www.uniprot.org/help/privacy).
my $agent = LWP::UserAgent->new( agent => "libwww-perl $contact" );
my $response = $agent->mirror( $query, $file );

if ( $response->is_success ) {
	my $results      = $response->header('X-Total-Results');
	my $release      = $response->header('X-UniProt-Release');
	my $release_date = $response->header('X-UniProt-Release-Date');
	print
"Downloaded proteome metadata for organism ID: $organism_id from UniProt release $release ($release_date) to file $file\n";
}
elsif ( $response->code == HTTP::Status::RC_NOT_MODIFIED ) {
print "Data for taxon $organism_id is up-to-date.\n";
}
else {
	die 'Failed, got '
	. $response->status_line . ' for '
	. $response->request->uri . "\n";
}
