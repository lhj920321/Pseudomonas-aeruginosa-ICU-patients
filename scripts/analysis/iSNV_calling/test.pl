use warnings;
use strict;


my @seq = (11,11,91,0);

my ($ra_1, $ra_2) = get_des_rank(\@seq);

print "@seq\n";
print "A G C T\n";
print "@$ra_1, @$ra_2\n";

sub get_des_rank{
	my @arr = @{$_[0]};
	my $DEBUG = 1;
	print "get_des_rank();@arr\n" if ($DEBUG == 1);
	my @bases = qw/A G C T/;
	my $rh_oldrank;
	for(my $i = 0; $i < @bases; $i ++){
		$rh_oldrank->{$bases[$i]} = $i;
	}
	my $rh_val2key;
	for(my $i = 0; $i < @arr; $i++){
		if (defined $rh_val2key->{"$arr[$i]"}){
			$rh_val2key->{"$arr[$i]"} .= "/".$bases[$i];
		} else {
			$rh_val2key->{"$arr[$i]"} = $bases[$i];
		}
	}
	my @vals = sort {$b<=>$a} keys %$rh_val2key;
	my @newrank = ();
	my @newrankAlle = ();

	print "get_des_rank():@vals\n" if ($DEBUG == 1);
	for (my $i = 0; $i < @vals; $i++ ) {
		my $alle = $rh_val2key->{$vals[$i]};
		push @newrankAlle, $alle;
		# push the index
		print "get_des_rank():$alle\n" if ($DEBUG == 1);
		if (length($alle) == 1){
			push @newrank,$rh_oldrank->{$alle};
		} elsif (length($alle) > 1) {
			my @bs = split /\//,$alle;
			push @newrank, $rh_oldrank->{$bs[0]};
			push @newrank, $rh_oldrank->{$bs[1]};
			# *****************************************************************************
			# *    下面想处理简并碱基的情况，但是没有写完，直接赋值bs[0]了，待修改        *
			print "get_des_rank():-->@bs\n" if ($DEBUG == 1);
			my $degranks = '';
			foreach my $ba (@bs){
				die "get_des_rank(): cannot find rank for base $ba\n" if (! defined $rh_oldrank->{$ba});
				$degranks .= '/'.$rh_oldrank->{$ba};
			}
			print "degrank = $degranks\n";
			$degranks =~ s/^\///;
			$degranks = $bs[0]; # <---------- for simple in its parent sub, just get the first
			#push @newrank, $rh_oldrank->{$degranks};
			# *****************************************************************************

		}
	} 
	print "get_des_rank():@newrank\n" if ($DEBUG == 1);
	print "get_des_rank():@newrankAlle\n" if ($DEBUG == 1);
	@newrank = @newrank[0..1];
	return (\@newrank, \@newrankAlle);
}
