use 5.008006;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME              => 'Bio::Homology::InterologWalk',
    VERSION_FROM      => 'lib/Bio/Homology/InterologWalk.pm', # finds $VERSION
    PREREQ_PM         => {
                Bio::Perl                  => 0,
                Bio::EnsEMBL::Registry     => 0,
		GO::Parser                 => 0,
		REST::Client               => 0,
		DBD::CSV                   => 0,
		String::Approx             => 0,
		List::Util                 => 0,
                File::Glob                 => 0,
                Carp                       => 0			

	},
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM  => 'lib/Bio/Homology/InterologWalk.pm', # retrieve abstract from module
       AUTHOR         => 'Giuseppe Gallone <ggallone@cpan.org>') : ()),
);
