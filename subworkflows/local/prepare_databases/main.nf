include { ARIA2 as TAXDUMP_DOWNLOAD                         } from '../../../modules/nf-core/aria2/main'
include { BLAST_BLASTDBCMD as BLASTDBCMD_EXTRACT_FASTA      } from '../../../modules/nf-core/blast/blastdbcmd/main'
include { BLAST_BLASTDBCMD as BLASTDBCMD_EXTRACT_TAXMAPPING } from '../../../modules/nf-core/blast/blastdbcmd/main'
include { BLAST_UPDATEBLASTDB                               } from '../../../modules/nf-core/blast/updateblastdb/main'
include { MMSEQS_CREATEDB                                   } from '../../../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_CREATEINDEX                                } from '../../../modules/nf-core/mmseqs/createindex/main'
include { MMSEQS_CREATETAXDB                                } from '../../../modules/nf-core/mmseqs/createtaxdb/main'
include { MMSEQS_DATABASES                                  } from '../../../modules/nf-core/mmseqs/databases/main'
include { OSFCLIENT_FETCH as EMU_DB_FETCH                   } from '../../../modules/nf-core/osfclient/fetch/main'
include { UNTAR as EMU_DB_UNTAR                             } from '../../../modules/nf-core/untar/main'
include { UNTAR as TAXDUMP_UNTAR                            } from '../../../modules/nf-core/untar/main'


workflow PREPARE_DATABASES {
    take:
    emu_db_name
    mmseqs2_db_name
    val_skip_emu
    val_skip_mmseqs2

    main:
    ch_emu_db = channel.empty()
    ch_mmseqs2_db = channel.empty()

    if (!val_skip_emu) {
        db_paths = [
            'default': 'osfstorage/emu-prebuilt/emu.tar.gz',
            'rdp': 'osfstorage/emu-prebuilt/rdp.tar.gz',
            'silva': 'osfstorage/emu-prebuilt/silva_database.tar.gz',
        ]

        EMU_DB_FETCH(
            [
                [id: "emu_database_${emu_db_name}"],
                '56uf7',
                db_paths[emu_db_name],
            ]
        )

        EMU_DB_UNTAR(EMU_DB_FETCH.out.download_files)

        ch_emu_db = EMU_DB_UNTAR.out.untar
    }

    if (!val_skip_mmseqs2) {
        if (mmseqs2_db_name == 'genbank') {
            TAXDUMP_DOWNLOAD(
                [
                    [id: 'taxdump'],
                    'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz',
                ]
            )

            TAXDUMP_UNTAR(TAXDUMP_DOWNLOAD.out.downloaded_file)

            BLAST_UPDATEBLASTDB([[id: '16S_ribosomal_RNA'], '16S_ribosomal_RNA'])

            BLASTDBCMD_EXTRACT_FASTA(
                [[id: '16S_genbank'], 'all', []],
                BLAST_UPDATEBLASTDB.out.db,
            )

            BLASTDBCMD_EXTRACT_TAXMAPPING(
                [[id: '16S_genbank.acc2taxid'], 'all', []],
                BLAST_UPDATEBLASTDB.out.db,
            )

            MMSEQS_CREATEDB(BLASTDBCMD_EXTRACT_FASTA.out.fasta)

            MMSEQS_CREATEINDEX(MMSEQS_CREATEDB.out.db)

            MMSEQS_CREATETAXDB(
                MMSEQS_CREATEINDEX.out.db_indexed,
                TAXDUMP_UNTAR.out.untar,
                BLASTDBCMD_EXTRACT_TAXMAPPING.out.text,
            )

            ch_mmseqs2_db = MMSEQS_CREATETAXDB.out.db_with_taxonomy
        }
        else if (mmseqs2_db_name == 'silva') {
            MMSEQS_DATABASES('silva')

            ch_mmseqs2_db = MMSEQS_DATABASES.out.databases.map { db -> [[id: '16S_SILVA'], db] }
        }
    }

    emit:
    emu_database     = ch_emu_db
    mmseqs2_database = ch_mmseqs2_db
}
