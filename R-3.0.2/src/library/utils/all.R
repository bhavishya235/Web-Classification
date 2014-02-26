#  File src/library/utils/R/MARC.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


MARC_relator_db <-
structure(list(term = c("Actor", "Adapter", "Annotator", "Applicant",
"Architect", "Arranger", "Artist", "Assignee", "Associated name",
"Attributed name", "Auctioneer", "Author", "Author in quotations or text extracts",
"Author of afterword, colophon, etc.", "Author of dialog", "Author of introduction, etc.",
"Author of screenplay, etc.", "Bibliographic antecedent", "Binder",
"Binding designer", "Book designer", "Book producer", "Bookjacket designer",
"Bookplate designer", "Bookseller", "Calligrapher", "Cartographer",
"Censor", "Choreographer", "Client", "Collaborator", "Collector",
"Collotyper", "Commentator", "Commentator for written text",
"Compiler", "Complainant", "Complainant-appellant", "Complainant-appellee",
"Composer", "Compositor", "Conceptor", "Conductor", "Consultant",
"Consultant to a project", "Contestant", "Contestant-appellant",
"Contestant-appellee", "Contestee", "Contestee-appellant", "Contestee-appellee",
"Contractor", "Contributor", "Copyright claimant", "Copyright holder",
"Corrector", "Correspondent", "Costume designer", "Creator",
"Curator of an exhibition", "Dancer", "Dedicatee", "Dedicator",
"Defendant", "Defendant-appellant", "Defendant-appellee", "Degree grantor",
"Delineator", "Depositor", "Designer", "Director", "Dissertant",
"Distributor", "Donor", "Draftsman", "Dubious author", "Editor",
"Electrotyper", "Engineer", "Engraver", "Etcher", "Expert", "Facsimilist",
"Film editor", "Forger", "Former owner", "Funder", "Honoree",
"Host", "Illuminator", "Illustrator", "Inscriber", "Instrumentalist",
"Interviewee", "Interviewer", "Inventor", "Landscape architect",
"Lender", "Libelant", "Libelant-appellant", "Libelant-appellee",
"Libelee", "Libelee-appellant", "Libelee-appellee", "Librettist",
"Licensee", "Licensor", "Lithographer", "Lyricist", "Metadata contact",
"Metal-engraver", "Moderator", "Monitor", "Musician", "Narrator",
"Opponent", "Organizer of meeting", "Originator", "Other", "Owner",
"Papermaker", "Patent applicant", "Patent holder", "Patron",
"Performer", "Photographer", "Plaintiff", "Plaintiff-appellant",
"Plaintiff-appellee", "Platemaker", "Printer", "Printer of plates",
"Process contact", "Producer", "Production personnel", "Programmer",
"Proofreader", "Publisher", "Publishing director", "Recipient",
"Recording engineer", "Redactor", "Renderer", "Research team head",
"Research team member", "Researcher", "Respondent", "Respondent-appellant",
"Respondent-appellee", "Reviewer", "Rubricator", "Scenarist",
"Scientific advisor", "Scribe", "Sculptor", "Secretary", "Signer",
"Singer", "Speaker", "Sponsor", "Standards body", "Stereotyper",
"Surveyor", "Thesis advisor", "Transcriber", "Translator", "Type designer",
"Typographer", "Vocalist", "Witness", "Wood-engraver", "Woodcutter",
"Writer of accompanying material"), code = c("act", "adp", "ann",
"app", "arc", "arr", "art", "asg", "asn", "att", "auc", "aut",
"aqt", "aft", "aud", "aui", "aus", "ant", "bnd", "bdd ", "bkd",
"bkp", "bjd", "bpd", "bsl", "cll", "ctg", "cns", "chr", "cli",
"clb", "col", "clt", "cmm", "cwt", "com", "cpl", "cpt", "cpe",
"cmp", "cmt ", "ccp", "cnd", "csl", "csp", "cos", "cot", "coe",
"cts", "ctt", "cte", "ctr", "ctb", "cpc", "cph", "crr", "crp",
"cst", "cre", "cur", "dnc", "dte", "dto", "dfd", "dft", "dfe",
"dgg", "dln", "dpt", "dsr", "drt", "dis", "dst", "dnr", "drm",
"dub", "edt", "elt", "eng", "egr", "etr", "exp", "fac", "flm",
"frg ", "fmo", "fnd", "hnr", "hst", "ilu", "ill", "ins", "itr",
"ive", "ivr", "inv", "lsa", "len", "lil", "lit", "lie", "lel",
"let", "lee", "lbt", "lse", "lso", "ltg", "lyr", "mdc", "mte",
"mod", "mon", "mus", "nrt", "opn", "orm", "org", "oth", "own",
"ppm", "pta", "pth", "pat", "prf", "pht", "ptf", "ptt", "pte",
"plt", "prt", "pop", "prc", "pro", "prd", "prg", "pfr", "pbl",
"pbd", "rcp", "rce", "red", "ren", "rth", "rtm", "res", "rsp",
"rst", "rse", "rev", "rbr", "sce", "sad", "scr", "scl", "sec",
"sgn", "sng", "spk", "spn", "stn", "str", "srv", "ths", "trc",
"trl", "tyd", "tyg", "voc", "wit", "wde", "wdc", "wam"), description = c("Use for a person who principally exhibits acting skills in a musical or dramatic presentation or entertainment.",
"Use for a person who 1) reworks a musical composition, usually for a different medium, or 2) rewrites novels or stories for motion pictures or other audiovisual medium.",
"Use for a person who writes manuscript annotations on a printed item.",
"Appraiser (USE: Expert)", "", "Use for a person who transcribes a musical composition, usually for a different medium from that of the original, in an arrangement the musical substance remains essentially unchanged.",
"Use for a person (e.g., a painter) who conceives, and perhaps also implements, an original graphic design or work of art, if specific codes (e.g., [egr], [etr]) are not desired. For book illustrators, prefer Illustrator [ill]. (UF: Graphic technician)",
"Use for a person or organization to whom a license for printing or publishing has been transferred.",
"Use as a general relator for a name associated with or found in an item or collection, or which cannot be determined to be that of a Former owner [fmo] or other designated relator indicative of provenance.",
"Use to relate an author, artist, etc. to a work for which there is or once was substantial authority for designating that person as author, creator, etc. of the work. (UF: Supposed name)",
"Use for a person or corporate body in change or the estimation and public auctioning of goods, particularly books, artistic works, etc.",
"Use for a person or corporate body chiefly responsible for the intellectual or artistic content of a work, usually printed text. This term may also be used when more than one person or body bears such responsibility. (UF: Joint author)",
"Use for a person whose work is largely quoted or extracted in a works to which he or she did not contribute directly. Such quotations are found particularly in exhibition catalogs, collections of photographs, etc.",
"Use for a person or corporate body responsible for an afterword, postface, colophon, etc. but who is not the chief author of a work.",
"Use for a person or corporate body responsible for the dialog or spoken commentary for a screenplay or sound recording.",
"Use for a person or corporate body responsible for an introduction, preface, foreword, or other critical introductory matter, but who is not the chief author.",
"Use for a person or corporate body responsible for a motion picture screenplay, dialog, spoken commentary, etc.",
"Use for the author responsible for a work upon which the work represented by the catalog record is based. This may be appropriate for adaptations, sequels, continuations, indexes, etc.",
"", "(UF: Designer of binding)", "Use for the person or firm responsible for the entire graphic design of a book, including arrangement of type and illustration, choice of materials, and process used. (UF: Designer of book)",
"Use for the person or firm responsible for the production of books and other print media, if specific codes (e.g., [bkd], [egr], [tyd], [prt]) are not desired. (UF: Producer of book)",
"(UF: Designer of bookjacket)", "(UF: Designer of bookplate)",
"Bowdlerizer (USE: Censor)", "", "", "Use for a censor, bowdlerizer, expurgator, etc., official or private. (UF: Bowdlerizer, Expurgator)",
"Use for a person who composes or arranges dances or other movements (e.g., 'master of swords') for a musical or dramatic presentation or entertainment.",
"Use for a person or organization for whom another person or organization is acting.",
"Use for a person or corporate body that takes a limited part in the elaboration of a work of another person or corporate body that brings complements (e.g., appendices, notes) to the work.",
"Use for a person who has brought together material from various sources, which has been arranged, described, and cataloged as a collection. The collector is neither the creator of the material nor the person to whom manuscripts in the collection may have been addressed.",
"", "Use for a person who provides interpretation, analysis, or a discussion of the subject matter on a recording, motion picture, or other audiovisual medium.",
"Use for a person or corporate body responsible for the commentary or explanatory notes about a text. For the writer of manuscript annotations in a printed book, use Annotator [ann].",
"Use for a person who produces a work or publication by selecting and putting together material from the works of various persons or bodies.",
"Use for the party who applies to the courts for redress, usually in an equity proceeding.",
"Use for a complainant who takes an appeal from one court or jurisdiction to another to reverse the judgment, usually in an equity proceeding.",
"Use for a complainant against whom an appeal is taken from one court or jurisdiction to another to reverse the judgment, usually in an equity proceeding.",
"Use for a person who creates a musical work, usually a piece of music in manuscript or printed form.",
"(UF: Typesetter)", "Use for a person or corporate body responsible for the original idea on which a work is based, this includes the scientific author of an audio-visual item and the conceptor of an advertisement.",
"Use for a person who directs a performing group (orchestra, chorus, opera, etc.).",
"Use for the person called upon for professional advice or services in a specialized field of knowledge or training.",
"Use for a person or corporate body engaged specifically to provide an intellectual overview of a strategic or operational task and by analysis, specification, or instruction, to create or propose a cost-effective course of action or solution.",
"Use for the party who opposes, resists, or disputes, in a court of law, a claim, decision, result, etc.",
"Use for a contestant who takes an appeal from one court of law or jurisdiction to another to reverse the judgment.",
"Use for a contestant against whom an appeal is taken from one court of law or jurisdiction to another to reverse the judgment.",
"Use for the party defending a claim, decision, result, etc. being opposed, resisted, or disputed in a court of law.",
"Use for a contestee who takes an appeal from one court or jurisdiction to another to reverse the judgment.",
"Use for a contestee against whom an appeal is taken from one court or jurisdiction to another to reverse the judgment.",
"Use for the person or corporate body who enters into a contract with another person or corporate body to perform a specific task.",
"Use for one whose work has been contributed to a larger work, such as an anthology, serial publication, or other compilation of individual works. Do not use for someone whose sole function in relation to a work is as author, editor, compiler or translator.",
"Use for the person listed as a copyright owner at the time of registration. Copyright can be granted or later transferred to another person or agent, at which time the claimant becomes the copyright holder.",
"", "Use for a corrector of manuscripts, such as the scriptorium official who corrected the work of a scribe. For printed matter, use Proofreader [pfr].",
"Use for a person or organization who was either the writer or recipient of a letter or other communication.",
"Use for a person who designs or makes costumes, fixes hair, etc., for a musical or dramatic presentation or entertainment. // Counterfeiter (USE: Forger)",
"Use for a person or corporate body responsible for the intellectual or artistic content of a work.",
"Use for a person who is responsible for conceiving and organizing an exhibition.",
"Use for a person who principally exhibits dancing skills in a musical or dramatic presentation or entertainment.",
"Use for a person or organization to whom a book, manuscript, etc., is dedicated (not the recipient of a gift).",
"Use for the author of a dedication, which may be a formal statement or in epistolary or verse form.",
"Use for the party defending or denying allegations made in a suit and against whom relief or recovery is sought in the courts, usually in a legal action.",
"Use for a defendant who takes an appeal from one court or jurisdiction to another to reverse the judgment, usually in a legal action.",
"Use for a defendant against whom an appeal is taken from one court or jurisdiction to another to reverse the judgment, usually in a legal action.",
"Use for the corporate body granting a degree for which the thesis or dissertation described was presented.",
"Use for a person or organization executing technical drawings from others' designs. // Deponent (USE: Witness)",
"Use for a person or organization placing material in the physical custody of a library or repository without transferring the legal title.",
"Use for a person or organization responsible for design if specific codes (e.g., [bkd], [tyd]) are not desired. // Designer of binding (USE: Binding designer) // Designer of book (USE: Book designer) // Designer of bookjacket (USE: Bookjacket designer) // Designer of bookplate (USE: Bookplate designer) // Designer of type (USE: Type designer)",
"Use for a person who is responsible for the general management of a work or who supervises the production of a performance for stage, screen, or sound recording.",
"Use for a person who presents a thesis for a university or higher-level educational degree.",
"Use for an agent or agency that has exclusive or shared marketing rights for an item.",
"Use for the donor of a book, manuscript, etc., to its present owner. Donors to previous owners are designated as Former owner [fmo] or Inscriber [ins].",
"Use for the person who prepares technical or mechanical drawings. (UF: Technical draftsman)",
"Use for a person or corporate body to which authorship has been dubiously or incorrectly ascribed.",
"Use for a person who prepares for publication a work not primarily his/her own, such as by elucidating text, adding introductory or other critical matter, or technically directing an editorial staff.",
"", "Use for a person or organization that is responsible for technical planning and design, particularly with construction.",
"", "", "Use for a person in charge of the description and appraisal of the value of goods, particularly rare items, works of art, etc. (UF: Appraiser) // Eyewitness (USE: Witness) // Expurgator (USE: Censor)",
"Use for the person or body that executed the facsimile. (UF: Copier)",
"Use for an editor of a motion picture film. This term is used regardless of the medium upon which the motion picture is produced or manufactured (e.g., acetate film, video tape). (UF: Motion picture editor)",
"(UF: Copier, Counterfeiter)", "Use for the person or organization who owned an item at any time in the past. Includes those to whom the material was once presented. The person or organization giving the item to the present owner is designated as Donor [dnr]",
"Use for the person or agency that furnished financial support for the production of the work. // Graphic technician (USE: Artist) [Relator term 'Graphic technician' (coded [grt]) used before March 1988 only.]",
"Use for the person in memory or honor of whom a book, manuscript, etc. is donated. (UF: Memorial)",
"Use for the person who is invited or regularly leads a program (often broadcast) that includes other guests, performers, etc. (e.g., talk show host).",
"", "Use for the person who conceives, and perhaps also implements, a design or illustration, usually to accompany a written text. // Imprimatur (USE: Licensor)",
"Use for the person who signs a presentation statement.", "Use for a person who principally plays an instrument in a musical or dramatic presentation or entertainment.",
"", "", "// Investigator (USE: Originator) // Joint author (USE: Author)",
"Use for the person or organization whose work involves coordinating the arrangement of existing and proposed land features and structures. ",
"Use for a person or organization permitting the temporary use of a book, manuscript, etc., such as for photocopying or microfilming.",
"Use for the party who files a libel in an ecclesiastical or admiralty case.",
"Use for a libelant who takes an appeal from one ecclesiastical court or admiralty to another to reverse the judgment.",
"Use for a libelant against whom an appeal is taken from one ecclesiastical court or admiralty to another to reverse the judgment.",
"Use for the party against whom a libel has been filed in an ecclesiastical court or admiralty.",
"Use for a libelee who takes an appeal from one ecclesiastical court or admiralty to another to reverse the judgment.",
"Use for a libelee against whom an appeal is taken from one ecclesiastical court or admiralty to another to reverse the judgment.",
"Use for the writer of the text of an opera, oratorio, etc.",
"Use for the original recipient of the right to print or publish.",
"Use for the signer of the license, imprimatur, etc. (UF: Imprimatur)",
"Use for the person who prepares the stone or plate for lithographic printing, including a graphic artist creating a design directly on the surface from which printing will be done.",
"Use for the writer of the text of a song. // Memorial (USE: Honoree)",
"Use for the person or organization primarily responsible for compiling and maintaining the original description of a metadata set (e.g., geospatial metadata set).",
"", "Use for the person who leads a program (often broadcast) where topics are discussed, usually with participation of experts in fields related to the discussion.",
"Use for a person or organization that supervises compliance with the contract and is responsible for the report and controls its distribution. Sometimes referred to as the grantee, or controlling agency. // Motion picture editor (USE: Film editor)",
"Use for the person who performs music or contributes to the musical content of a work when it is not possible or desirable to identify the function more precisely.",
"Use for the speaker who relates the particulars of an act, occurrence, or course of events. // Observer (USE: Witness) // Onlooker (USE: Witness)",
"Use for the person or corporate body responsible for opposing a thesis or dissertation.",
"Use for the person or corporate body responsible for organizing a meeting for which an item is the report or proceedings.",
"Use for the author or agency performing the work, i.e., the name of a person or organization associated with the intellectual content of the work. This category does not include the publisher or personal affiliation, or sponsor except where it is also the corporate author. Includes a person designated in the work as investigator or principal investigator. (UF: Principal investigator)",
"Use for relator codes from other lists which have no equivalent in the MARC list or for terms which have not been assigned a code.",
"Use for the person or organization that currently owns an item or collection.",
"", "Use for the person or corporate body that applied for a patent.",
"Use for the person or corporate body that was granted the patent referred to by the item. (UF: Patentee) // Patentee (USE: Patent holder)",
"Use for the person responsible for commissioning a work. Usually a patron uses his or her means or influence to support the work of artists, writers, etc. This includes those who commission and pay for individual works.",
"User for a person who exhibits musical or acting skills i a musical or dramatic presentation or entertainment, if specific codes for those functions ([act], [dnc], [itr], [voc], etc.) are not used. If specific codes are used, [prf] is used for a person whose principal skill is not known or specified. // Performer of research (USE: Researcher)",
"Use for the person or organization responsible for taking photographs, whether they are used in their original form or as reproductions.",
"Use for the party who complains or sues in court in a personal action, usually in a legal proceeding.",
"Use for a plaintiff who takes an appeal from one court or jurisdiction to another to reverse the judgment, usually in a legal proceeding.",
"Use for a plaintiff against whom an appeal is taken from one court or jurisdiction to another to reverse the judgment, usually in a legal proceeding.",
"// Plates, Printer of (USE: Printer of Plates) // Principal investigator (USE: Originator)",
"Use for the person or organization who prints texts, whether from type or plates.",
"Use for the person or organization who prints illustrations from plates. (UF: Plates, Printer of)",
"Use for a person or organization primarily responsible for performing or initiating a process, such as is done with the collection of metadata sets.",
"Use for a person who is responsible for the making of a motion picture, including business aspects, management of the productions, and the commercial success of the work. // Producer of book (USE: Book producer)",
"Use for a person who is associated with the production (props, lighting, special effects, etc.) of a musical or dramatic presentation or entertainment.",
"Use for a person or corporate body responsible for the creation and/or maintenance of computer program design documents, source code, and machine-executable digital files and supporting documentation. // Promoter (USE: Thesis advisor)",
"Use for a person who corrects printed matter. For manuscripts, use Corrector [crr].",
"", "Use for a person who presides over the elaboration of a collective work to ensure its coherence or continuity. This includes editors-in-chief, literary editors, editors of series, etc.",
"Use for the person to whom correspondence is addressed.", "Use for a person who supervises the technical aspects of a sound or video recording session.",
"Use for a person who writes or develops the framework for an item without being intellectually responsible for its content.",
"Use for the draftsman who prepares drawings of architectural designs (i.e., renderings) in accurate, representational perspective to show what the project will look like when completed.",
"Use for the person or corporate body that directed or managed a research project.",
"Use for the person or corporate body that participated in a research project but whose role did not involve direction or management of it.",
"Use for the person or corporate body responsible for performing research. (UF: Performer of research)",
"Use for the party who makes an answer to the courts pursuant to an application for redress, usually in an equity proceeding.",
"Use for a respondent who takes an appeal from one court or jurisdiction to another to reverse the judgment, usually in an equity proceeding.",
"Use for a respondent against whom an appeal is taken from one court or jurisdiction to another to reverse the judgment, usually in an equity proceeding.",
"Use for a person or corporate body responsible for the review of book, motion picture, performance, etc.",
"", "Use for the author of a motion picture screenplay.", "Use for a person who brings scientific, pedagogical, or historical competence to the conception and realization on a work, particularly in the case of audio-visual items.",
"Use for an amanuensis and for a writer of manuscripts proper. For a person who makes pen-facsimiles, use Facsimilist [fac].",
"Use when the more general term Artist [art] is not desired.",
"Use for a recorder, redactor, or other person responsible for expressing the views of a corporate body.",
"Use for the person whose signature appears without a presentation or other statement indicative of provenance. When there is a presentation statement, use Inscriber [ins].",
"Use for a person who uses his or her voice with or without instrumental accompaniment to produce music. A singer's performance may or may not include actual words.",
"Use for a person who participates in a program (often broadcast) and makes a formalized contribution or presentation generally prepared in advance.",
"Use for the person or agency that issued a contract or under the auspices of which a work has been written, printed, published, etc.",
"Use for a corporate body or agency responsible for the development or enforcement of a standard.",
"// Supposed name (USE: Attributed name)", "Use for a person or organization who does measurements of tracts of land, etc. to determine location, forms, and boundaries. // Technical draftsman (USE: Draftsman) // Testifier (USE: Witness)",
"Use for the person under whose supervision a degree candidate develops and presents a thesis, memoire, or text of a dissertation. (UF: Promoter)",
"Use for a person who prepares a handwritten or typewritten copy from original material, including from dictated or orally recorded material. For makers of pen-facsimiles, use Facsimilist [fac].",
"Use for a person who renders a text from one language into another, or from an older form of a language into the modern form.",
"Use for the person who designed the type face used in a particular item. (UF: Designer of type) // Typesetter (USE: Compositor)",
"Use for the person primarily responsible for choice and arrangement of type used in an item. If the typographer is also responsible for other aspects of the graphic design of a book (e.g., Book designer [bkd]), codes for both functions may be needed.",
"Use for a person who principally exhibits singing skills in a musical or dramatic presentation or entertainment.",
"Use for a person who verifies the truthfulness of an event or action. (UF: Deponent, Eyewitness, Observer, Onlooker, Testifier)",
"User for a person who makes prints by cutting the image in relief on the end-grain of a wood block.",
"User for a person who makes prints by cutting the image in relief on the plank side of a wood block.",
"Use for a person who writes significant material which accompanies a sound recording or other audiovisual material."
), usage = c("", "", "", "", "", "", "", "", "", "", "", "Use for full authors who have made substantial contributions to the package and should show up in the package citation.",
"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
"", "", "", "", "", "", "", "Use for package maintainers that collected code (potentially in other languages) but did not make further substantial contributions to the package.",
"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
"Use for authors who have made smaller contributions (such as code patches etc.) but should not show up in the package citation.",
"", "Use for all copyright holders.", "", "", "", "Use for the package maintainer.",
"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
"", "", "", "", "", "", "", "", "If the package is part of a thesis, use for the thesis advisor.",
"", "If the R code is merely a translation from another language (typically S), use for the translator to R.",
"", "", "", "", "", "", "")), .Names = c("term", "code", "description",
"usage"), class = "data.frame", row.names = c(NA, -173L))
MARC_relator_db_codes_used_with_R <-
c("aut", "com", "ctb", "cph", "cre", "ths", "trl")
#  File src/library/utils/R/RShowDoc.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

RShowDoc <- function(what, type=c("pdf", "html", "txt"), package)
{
    paste. <- function(x, ext) paste(x, ext, sep=".")
    pdf_viewer <- function(path) {
        pdfviewer <- getOption("pdfviewer")
        if(identical(pdfviewer, "false")) {
        } else if(.Platform$OS.type == "windows" &&
                  identical(pdfviewer, file.path(R.home("bin"), "open.exe")))
            shell.exec(path)
        else system2(pdfviewer, shQuote(path), wait = FALSE)
    }

    html_viewer <- function(path) {
        ## we don't use browseURL under Windows as shell.exec does
        ## not want an encoded URL.
        browser <- getOption("browser")
        if(is.null(browser) && .Platform$OS.type == "windows")
            shell.exec(chartr("/", "\\", path))
        else browseURL(paste0("file://", URLencode(path)))
    }

    type <- match.arg(type)
    if(missing(what) || length(what) != 1L || !is.character(what)) {
        message("   RShowDoc() should be used with a character string argument specifying\n   a documentation file")
        return(invisible())
    }
    if(!missing(package)) {
        pkgpath <- find.package(package)
        if(type == "pdf") {
            path <- file.path(pkgpath, "doc", paste.(what, "pdf"))
            if(file.exists(path)) {
                pdf_viewer(path)
                return(invisible(path))
            }
            path <- file.path(pkgpath, paste.(what, "pdf"))
            if(file.exists(path)) {
                pdf_viewer(path)
                return(invisible(path))
            }
            type <- "html"
        }
        if(type == "html") {
            path <- file.path(pkgpath, "doc", paste.(what, "html"))
            if(file.exists(path)) {
                html_viewer(path)
                return(invisible(path))
            }
            path <- file.path(pkgpath, paste.(what, "html"))
            if(file.exists(path)) {
                html_viewer(path)
                return(invisible(path))
            }
        }
        path <- file.path(pkgpath, "doc", what)
        if(file.exists(path)) {
            file.show(path)
            return(invisible(path))
        }
        path <- file.path(pkgpath, what)
        if(file.exists(path)) {
            file.show(path)
            return(invisible(path))
        }
        stop(gettextf("no documentation for %s found in package %s",
                      sQuote(what), sQuote(package)), domain = NA)
    }
    if(what == "FAQ") what <- "R-FAQ"
    if(what == "NEWS") {
        if(type == "pdf") type <- "html"
        if(type == "html") {
            path <- file.path(R.home("doc"), "html", paste.(what, "html"))
            if(file.exists(path)) {
                html_viewer(path)
                return(invisible(path))
            }
        }
        ## This is in UTF-8 and has a BOM on the first line
        path <- file.path(R.home(), what)
        tf <- tempfile()
        tmp <- readLines(path)
        tmp[1] <- ""
        writeLines(tmp, tf)
        file.show(tf, delete.file = TRUE, encoding = "UTF-8")
        return(invisible(path))
    } else if(what == "COPYING") {
        path <- file.path(R.home("doc"), what)
        file.show(path)
        return(invisible(path))
    } else if(what %in% dir(file.path(R.home("share"), "licenses"))) {
        path <- file.path(R.home("share"), "licenses", what)
        file.show(path)
        return(invisible(path))
    } else if(what %in% c("R-admin", "R-data", "R-exts", "R-FAQ", "R-intro",
                          "R-ints", "R-lang")) {
        if(type == "pdf") {
            path <- file.path(R.home("doc"), "manual", paste.(what, "pdf"))
            if(file.exists(path)) {
                pdf_viewer(path)
                return(invisible(path))
            }
            type <- "html"
        }
        if(type == "html") {
            path <- file.path(R.home("doc"), "manual", paste.(what, "html"))
            if(file.exists(path)) {
                html_viewer(path)
                return(invisible(path))
            }
        }
        if(what == "R-FAQ" &&
           file.exists(path <- file.path(R.home("doc"), "FAQ"))) {
            file.show(path)
            return(invisible(path))
        }
    } else if(.Platform$OS.type == "windows" && what %in% "rw-FAQ") {
        if(type == "pdf") type <- "html"
        if(type == "html") {
            path <- file.path(R.home("doc"), "html", paste.(what, "html"))
            if(file.exists(path)) {
                html_viewer(path)
                return(invisible(path))
            }
        }
        path <- file.path(R.home("doc"), what)
        if(file.exists(path)) {
            file.show(path)
            return(invisible(path))
        }
        path <- file.path(R.home(), "src", "gnuwin32", what)
        if(file.exists(path)) {
            file.show(path)
            return(invisible(path))
        }
    } else {
        rdocdir <- R.home("doc")
        docs <- dir(rdocdir, full.names=TRUE)
        docs <- docs[sapply(docs, function(x) file_test("-f", x))]
        m <- match(what, basename(docs), 0L)
        if(m > 0L) {
            file.show(docs[m])
            return(invisible(docs[m]))
        }
    }
    stop("document not found")
}
#  File src/library/utils/R/RSiteSearch.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

RSiteSearch <- function(string,
                        restrict = c("functions", "vignettes", "views"),
			format = c("normal", "short"),
			sortby = c("score", "date:late", "date:early",
			"subject", "subject:descending",
			"from", "from:descending", "size", "size:descending"),
			matchesPerPage = 20)
{
    string <- paste0("http://search.r-project.org/cgi-bin/namazu.cgi?query=",
		     gsub(" ", "+", string))
    mpp <- paste0("max=", matchesPerPage)
    format <- paste0("result=", match.arg(format))

    restrictVALS <- c("functions", "vignettes", "views")
    restr <- match.arg(restrict, choices = restrictVALS, several.ok = TRUE)
    restr <- paste(paste0("idxname=", restr), collapse = "&")

    sortby <- match.arg(sortby)
    sortby <- paste0("sort=",
		     switch(sortby,
			    "score"=, "date:late"=, "date:early" = sortby,
			    "subject"		 = "field:subject:ascending",
			    "subject:descending" = "field:subject:descending",
			    "from"		 = "field:from:ascending",
			    "from:descending"	 = "field:from:descending",
			    "size"		 = "field:size:ascending",
			    "size:descending"	 = "field:size:descending"))

    ## we know this is a http:// URL, so encoding should be safe.
    ## it seems that firefox on Mac OS needs it for {...}
    ## OTOH, Namazu does not decode in, say, sort=date:late.
    qstring <- paste(URLencode(string), mpp, format, sortby, restr, sep = "&")
    browseURL(qstring)
    cat(gettextf("A search query has been submitted to %s",
                 "http://search.r-project.org"), "\n", sep = "")
    cat(gettext("The results page should open in your browser shortly\n"))
    invisible(qstring)
}
#  File src/library/utils/R/Rprof.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

Rprof <- function(filename = "Rprof.out", append = FALSE, interval =  0.02,
                  memory.profiling = FALSE, gc.profiling = FALSE,
                  line.profiling = FALSE, numfiles = 100L, bufsize = 10000L)
{
    if(is.null(filename)) filename <- ""
    invisible(.External(C_Rprof, filename, append, interval, memory.profiling,
                        gc.profiling, line.profiling, numfiles, bufsize))
}

Rprofmem <- function(filename = "Rprofmem.out", append = FALSE, threshold = 0)
{
    if(is.null(filename)) filename <- ""
    invisible(.External(C_Rprofmem, filename, append, as.double(threshold)))
}
#   File src/library/utils/R/Sweave.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

### The drivers are now in SweaveDrivers.R

### FIXMEs
### b) It would be nice to allow multiple 'grdevice' options

### Encodings (currently, different from 2.13.0)
###
### SweaveReadFile figures out an encoding, uses it (not currently for
### \SweaveInclude files) and returns it as an attribute.  This is
### then passed as an attribute of 'file' to the driver's setup
### routine.  Unless it is "" or "ASCII", the RweaveLatex driver
### re-encodes the output back to 'encoding': the Rtangle driver
### leaves it in the encoding of the current locale and records what
### that is in a comment.
###
### SweaveReadFile first looks for a call to one of the LaTeX packages
### inputen[cx] and deduces the vignette encoding from that, falling
### back to the package encoding, then Latin-1 (with a warning).  This
### should work OK provided the package encoding is Latin-1: it is
### UTF-8 then LaTeX needs to be told what to do.  It also assumes
### that R output is in the current locale: a package with a different
### encoding from the current one might have data in that package's
### encoding.


### Correspondence between input and output is maintained in two
### places: Each chunk has a srclines attribute, recording the input
### lines it corresponds to.  Each code chunk will have attached
### srcrefs that duplicate the srclines.  We don't need srclines for
### code, but we do need it for doc chunks, and it's easiest to just
### keep it for everything.



Stangle <- function(file, driver = Rtangle(),
                    syntax = getOption("SweaveSyntax"),
                    encoding = "", ...)
    Sweave(file = file, driver = driver, encoding = encoding, ...)

Sweave <- function(file, driver = RweaveLatex(),
                   syntax = getOption("SweaveSyntax"),
                   encoding = "", ...)
{
    if (is.character(driver)) driver <- get(driver, mode = "function")()
    else if (is.function(driver)) driver <- driver()

    if (is.null(syntax)) syntax <- SweaveGetSyntax(file) # from the extension
    if (is.character(syntax)) syntax <- get(syntax, mode = "list")

    if (.Platform$OS.type == "windows") file <- chartr("\\", "/", file)

    text <- SweaveReadFile(file, syntax, encoding = encoding)
    attr(file, "encoding") <- encoding <- attr(text, "encoding")
    srcFilenames <- attr(text, "files")
    srcFilenum <- attr(text, "srcFilenum")
    srcLinenum <- attr(text, "srcLinenum")

    ## drobj$options is the current set of options for this file.
    drobj <- driver$setup(file = file, syntax = syntax, ...)
    on.exit(driver$finish(drobj, error = TRUE))

    syntax <- attr(text, "syntax") # this is from the file commands.

    if (!is.na(envopts <- Sys.getenv("SWEAVE_OPTIONS", NA)))
        drobj$options <-
            SweaveParseOptions(envopts, drobj$options, driver$checkopts)

    drobj$filename <- file

    mode <- "doc"
    chunknr <- 0L
    chunk <- NULL
    chunkopts <- NULL

    namedchunks <- list()
    prevfilenum <- 0L
    prevlinediff <- 0L
    for (linenum in seq_along(text)) {
    	line <- text[linenum]
    	filenum <- srcFilenum[linenum]
    	linediff <- srcLinenum[linenum] - linenum
	if(nzchar(Sys.getenv("R_DEBUG_Sweave"))) {
	    ## Extensive logging for debugging, needs 'ls' (unix-like or Rtools):
	    cat(sprintf("l.%3d: %30s -'%4s'- ", linenum, substr(line,1,30), mode))
	    cat(sprintf("%16s\n", system(paste("ls -s",
				   summary(drobj$output)$description), intern=TRUE)))
	}
        if (length(grep(syntax$doc, line))) { # start new documentation chunk
            if (mode == "doc") {
                if (!is.null(chunk)) drobj <- driver$writedoc(drobj, chunk)
            } else {
                if (!is.null(chunkopts$label))
                    namedchunks[[chunkopts$label]] <- chunk
                if (!is.null(chunk))
                    drobj <- driver$runcode(drobj, chunk, chunkopts)
                mode <- "doc"
            }
            chunk <- NULL
        } else if (length(grep(syntax$code, line))) { # start new code chunk
            if (mode == "doc") {
                if (!is.null(chunk)) drobj <- driver$writedoc(drobj, chunk)
            } else {
                if (!is.null(chunkopts$label))
                    namedchunks[[chunkopts$label]] <- chunk
                if (!is.null(chunk))
                    drobj <- driver$runcode(drobj, chunk, chunkopts)
            }
            mode <- "code"
            chunkopts <- sub(syntax$code, "\\1", line)
            chunkopts <- SweaveParseOptions(chunkopts,
                                            drobj$options,
                                            driver$checkopts)
            ## these #line directives are used for error messages when parsing
            file <- srcFilenames[filenum]
            chunk <- paste0("#line ", linenum+linediff+1L, ' "', basename(file), '"')
            attr(chunk, "srclines") <- linenum + linediff
            attr(chunk, "srcFilenum") <- filenum
            attr(chunk, "srcFilenames") <- srcFilenames
            chunknr <- chunknr + 1L  # this is really 'code chunk number'
            chunkopts$chunknr <- chunknr
        } else {  # continuation of current chunk
            if (mode == "code" && length(grep(syntax$coderef, line))) {
                chunkref <- sub(syntax$coderef, "\\1", line)
                if (!(chunkref %in% names(namedchunks))) {
                    ## omit unknown references
                    warning(gettextf("reference to unknown chunk %s",
                                     sQuote(chunkref)),
                            call. = TRUE,domain = NA)
                    next
                } else {
                    ## these #line directives are used for error messages
                    ## when parsing
                    file <- srcFilenames[filenum]
                    line <- c(namedchunks[[chunkref]],
			      paste0("#line ", linenum+linediff+1L,
				     ' "', basename(file), '"'))
                }
            }
            if (mode == "code" &&
                (prevfilenum != filenum ||
                 prevlinediff != linediff)) {
                file <- srcFilenames[filenum]
                line <- c(paste0("#line ", linenum+linediff, ' "', basename(file), '"'),
                          line)
            }
            srclines <- c(attr(chunk, "srclines"), rep(linenum+linediff, length(line)))
            srcfilenum <- c(attr(chunk, "srcFilenum"), rep(filenum, length(line)))
	    chunk <- c(chunk, line)
            attr(chunk, "srclines") <- srclines
            attr(chunk, "srcFilenum") <- srcfilenum
            attr(chunk, "srcFilenames") <- srcFilenames
	}
	prevfilenum <- filenum
	prevlinediff <- linediff
    }
    if (!is.null(chunk)) { # write out final chunk
	drobj <-
	    if (mode == "doc") driver$writedoc(drobj, chunk)
	    else driver$runcode(drobj, chunk, chunkopts)
    }

    on.exit() # clear action to finish with error = TRUE
    drobj$srcFilenames <- srcFilenames
    driver$finish(drobj)
}

SweaveReadFile <- function(file, syntax, encoding = "")
{
    ## file can be a vector to keep track of recursive calls to
    ## SweaveReadFile.  In this case only the first element is
    ## tried to read in, the rest are forbidden names for further
    ## SweaveInput
    f <- file[1L]

    bf <- basename(f)
    df <- dirname(f)
    if (!file.exists(f)) {
        f <- list.files(df, full.names = TRUE,
                        pattern = paste0(bf, syntax$extension))

        if (length(f) == 0L)
            stop(gettextf("no Sweave file with name %s found",
                          sQuote(file[1L])), domain = NA)
        else if (length(f) > 1L)
            stop(paste(gettextf("%d Sweave files for basename %s found",
                                length(f), sQuote(file[1L])),
                       paste(":\n         ", f, collapse="")),
                 domain = NA)
    }

    ## An incomplete last line is not a real problem.
    text <- readLines(f[1L], warn = FALSE)
    srcLinenum <- seq_along(text)

    if (encoding != "bytes")  {
        ## now sort out an encoding, if needed.
        enc <- tools:::.getVignetteEncoding(text, convert = TRUE)
        if (enc == "non-ASCII") {
            enc <- if (nzchar(encoding)) {
                encoding
            } else {
                stop(sQuote(basename(file)),
                        " is not ASCII and does not declare an encoding",
                        domain = NA, call. = FALSE)
            }
        } else if (enc == "unknown") {
            stop(sQuote(basename(file)),
                 " declares an encoding that Sweave does not know about",
                 domain = NA, call. = FALSE)
        }
        if (nzchar(enc)) text <- iconv(text, enc, "") else enc <- "ASCII"
    } else enc <- "bytes"

    pos <- grep(syntax$syntaxname, text)

    if (length(pos) > 1L)
        warning(gettextf("more than one syntax specification found, using the first one"), domain = NA)

    if (length(pos) > 0L) {
        sname <- sub(syntax$syntaxname, "\\1", text[pos[1L]])
        syntax <- get(sname, mode = "list")
        if (!identical(class(syntax), "SweaveSyntax"))
            stop(gettextf("object %s does not have class \"SweaveSyntax\"",
                          sQuote(sname)), domain = NA)
        text <- text[-pos]
        srcLinenum <- srcLinenum[-pos]
    }
    srcFilenum <- rep(1, length(srcLinenum))

    if (!is.null(syntax$input)) {
        while(length(pos <- grep(syntax$input, text))) {
            pos <- pos[1L]
            ifile <- file.path(df, sub(syntax$input, "\\1", text[pos]))
            if (any(ifile == file)) {
                stop(paste(gettextf("recursive Sweave input %s in stack",
                                    sQuote(ifile)),
                           paste("\n         ", seq_len(file), ": ",
                                 rev(file), collapse="")),
                 domain = NA)
            }
            itext <- SweaveReadFile(c(ifile, file), syntax, encoding = encoding)

	    pre <- seq_len(pos-1L)
	    post <- seq_len(length(text) - pos) + pos
	    text <- c(text[pre], itext, text[post])

	    srcLinenum <- c(srcLinenum[pre], attr(itext, "srcLinenum"),
	    		    srcLinenum[post])
	    srcFilenum <- c(srcFilenum[pre], attr(itext, "srcFilenum")+length(f),
	    		    srcFilenum[post])
	    f <- c(f, attr(itext, "files"))
        }
    }

    attr(text, "syntax") <- syntax
    attr(text, "files") <- f
    attr(text, "encoding") <- enc
    attr(text, "srcLinenum") <- srcLinenum
    attr(text, "srcFilenum") <- srcFilenum
    text
}



###**********************************************************

SweaveSyntaxNoweb <-
    list(doc = "^@",
         code = "^<<(.*)>>=.*",
         coderef = "^<<(.*)>>.*",
         docopt = "^[[:space:]]*\\\\SweaveOpts\\{([^\\}]*)\\}",
         docexpr = "\\\\Sexpr\\{([^\\}]*)\\}",
         extension = "\\.[rsRS]?nw$",
         syntaxname = "^[[:space:]]*\\\\SweaveSyntax\\{([^\\}]*)\\}",
         input = "^[[:space:]]*\\\\SweaveInput\\{([^\\}]*)\\}",
         trans = list(
             doc = "@",
             code = "<<\\1>>=",
             coderef = "<<\\1>>",
             docopt = "\\\\SweaveOpts{\\1}",
             docexpr = "\\\\Sexpr{\\1}",
             extension = ".Snw",
             syntaxname = "\\\\SweaveSyntax{SweaveSyntaxNoweb}",
             input = "\\\\SweaveInput{\\1}")
         )

class(SweaveSyntaxNoweb) <- "SweaveSyntax"

SweaveSyntaxLatex <- SweaveSyntaxNoweb
SweaveSyntaxLatex$doc <-  "^[[:space:]]*\\\\end\\{Scode\\}"
SweaveSyntaxLatex$code <- "^[[:space:]]*\\\\begin\\{Scode\\}\\{?([^\\}]*)\\}?.*"
SweaveSyntaxLatex$coderef <- "^[[:space:]]*\\\\Scoderef\\{([^\\}]*)\\}.*"
SweaveSyntaxLatex$extension <- "\\.[rsRS]tex$"

SweaveSyntaxLatex$trans$doc <-  "\\\\end{Scode}"
SweaveSyntaxLatex$trans$code <- "\\\\begin{Scode}{\\1}"
SweaveSyntaxLatex$trans$coderef <- "\\\\Scoderef{\\1}"
SweaveSyntaxLatex$trans$syntaxname <- "\\\\SweaveSyntax{SweaveSyntaxLatex}"
SweaveSyntaxLatex$trans$extension <- ".Stex"

SweaveGetSyntax <- function(file)
{
    synt <- apropos("SweaveSyntax", mode = "list")
    for (sname in synt) {
        s <- get(sname, mode = "list")
        if (!identical(class(s), "SweaveSyntax")) next
        if (length(grep(s$extension, file))) return(s)
    }
    SweaveSyntaxNoweb
}


SweaveSyntConv <- function(file, syntax, output=NULL)
{
    if (is.character(syntax)) syntax <- get(syntax)

    if (!identical(class(syntax), "SweaveSyntax"))
        stop(gettextf("target syntax not of class %s",
                      dQuote("SweaveSyntax")),
             domain = NA)
    if (is.null(syntax$trans))
        stop("target syntax contains no translation table")

    insynt <- SweaveGetSyntax(file)
    text <- readLines(file)
    if (is.null(output))
        output <- sub(insynt$extension, syntax$trans$extension, basename(file))

    TN <- names(syntax$trans)

    for (n in TN)
        if (n != "extension") text <- gsub(insynt[[n]], syntax$trans[[n]], text)

    cat(text, file = output, sep = "\n")
    cat("Wrote file", output, "\n")
}


###**********************************************************

## parses an option string, from
## - the header of a code chunk
## - an \SweaveOpts{} statement (strangely, left to the drivers)
## - the value of environment variable SWEAVE_OPTIONS
##
## The format is name=value pairs with whitespace being discarded
## (and could have been done all at once).
SweaveParseOptions <- function(text, defaults = list(), check = NULL)
{
    x <- sub("^[[:space:]]*(.*)", "\\1", text)
    x <- sub("(.*[^[:space:]])[[:space:]]*$", "\\1", x)
    x <- unlist(strsplit(x, "[[:space:]]*,[[:space:]]*"))
    x <- strsplit(x, "[[:space:]]*=[[:space:]]*")

    ## only the first option may have no name: the chunk label
    if (length(x)) {
        if (length(x[[1L]]) == 1L) x[[1L]] <- c("label", x[[1L]])
    } else return(defaults)

    if (any(sapply(x, length) != 2L))
        stop(gettextf("parse error or empty option in\n%s", text), domain = NA)

    options <- defaults
    for (k in seq_along(x)) options[[ x[[k]][1L] ]] <- x[[k]][2L]

    ## This is undocumented
    if (!is.null(options[["label"]]) && !is.null(options[["engine"]]))
        options[["label"]] <-
            sub(paste0("\\.", options[["engine"]], "$"),
                "", options[["label"]])

    if (!is.null(check)) check(options) else options
}

## really part of the RweaveLatex and Rtangle drivers
SweaveHooks <- function(options, run = FALSE, envir = .GlobalEnv)
{
    if (is.null(SweaveHooks <- getOption("SweaveHooks"))) return(NULL)

    z <- character()
    for (k in names(SweaveHooks))
        if (nzchar(k) && is.logical(options[[k]]) && options[[k]])
            if (is.function(SweaveHooks[[k]])) {
                z <- c(z, k)
                if (run) eval(SweaveHooks[[k]](), envir=envir)
            }
    z # a character vector.
}

### For R CMD xxxx ------------------------------------------
.Sweave <- function(args = NULL)
{
    options(warn = 1)
    if (is.null(args)) {
        args <- commandArgs(TRUE)
        args <- paste(args, collapse=" ")
        args <- strsplit(args,'nextArg', fixed = TRUE)[[1L]][-1L]
    }

    Usage <- function() {
        cat("Usage: R CMD Sweave [options] file",
            "",
            "A front-end for Sweave",
            "",
            "Options:",
            "  -h, --help      print this help message and exit",
            "  -v, --version   print version info and exit",
            "  --driver=name   use named Sweave driver",
            "  --encoding=enc  default encoding 'enc' for file",
            "  --options=      comma-separated list of Sweave options",
            "  --pdf           convert to PDF document",
            "  --compact=      try to compact PDF document:",
            '                  "no" (default), "qpdf", "gs", "gs+qpdf", "both"',
            "  --compact       same as --compact=qpdf",
            "",
            "Report bugs at bugs.r-project.org .",
            sep = "\n")
    }
    do_exit <- function(status = 0L)
        q("no", status = status, runLast = FALSE)

    if (!length(args)) {
        Usage()
        do_exit(1L)
    }
    file <- character()
    driver <- encoding <- options <- ""
    toPDF <- FALSE
    compact <- Sys.getenv("_R_SWEAVE_COMPACT_PDF_", "no")
    while(length(args)) {
        a <- args[1L]
        if (a %in% c("-h", "--help")) {
            Usage()
            do_exit()
        }
        else if (a %in% c("-v", "--version")) {
            cat("Sweave front-end: ",
                R.version[["major"]], ".",  R.version[["minor"]],
                " (r", R.version[["svn rev"]], ")\n", sep = "")
            cat("",
                "Copyright (C) 2006-2011 The R Core Team.",
                "This is free software; see the GNU General Public License version 2",
                "or later for copying conditions.  There is NO warranty.",
                sep = "\n")
            do_exit()
        } else if (substr(a, 1, 9) == "--driver=") {
            driver <- substr(a, 10, 1000)
        } else if (substr(a, 1, 11) == "--encoding=") {
            encoding <- substr(a, 12, 1000)
        } else if (substr(a, 1, 10) == "--options=") {
            options <- substr(a, 11, 1000)
        } else if (a == "--pdf") {
            toPDF <- TRUE
        } else if (substr(a, 1, 10) == "--compact=") {
            compact <- substr(a, 11, 1000)
        } else if (a == "--compact") {
            compact <- "qpdf"
        } else if (substr(a, 1, 1) == "-") {
            message(gettextf("Warning: unknown option %s", sQuote(a)),
                    domain = NA)
        } else file <- c(file, a)
       args <- args[-1L]
    }
    if(length(file) != 1L) {
        Usage()
        do_exit(1L)
    }
    args <- list(file)
    if(nzchar(driver)) args <- c(args, driver)
    args <- c(args, encoding = encoding)
    if(nzchar(options)) {
        opts <- eval(parse(text = paste("list(", options, ")")))
        args <- c(args, opts)
    }
    do.call(Sweave, args)
    if (toPDF) {
        texfile <- basename(sub("\\.[rsRS][[:alpha:]]+$", ".tex", file))
        tools::texi2pdf(texfile, clean = TRUE)
        ofile <- sub("\\.tex$", ".pdf", texfile)
        message(gettextf("Created PDF document %s", sQuote(ofile)),
                domain = NA)
        if(compact != "no") {
            ## <NOTE>
            ## Same code as used for --compact-vignettes in
            ## .build_packages() ...
            message("Compacting PDF document")
            if(compact %in% c("gs", "gs+qpdf", "both")) {
                gs_cmd <- tools:::find_gs_cmd(Sys.getenv("R_GSCMD", ""))
                gs_quality <- "ebook"
            } else {
                gs_cmd <- ""
                gs_quality <- "none"
            }
            qpdf <- if(compact %in% c("qpdf", "gs+qpdf", "both"))
                Sys.which(Sys.getenv("R_QPDF", "qpdf"))
            else ""
            res <- tools::compactPDF(ofile, qpdf = qpdf,
                                     gs_cmd = gs_cmd,
                                     gs_quality = gs_quality)
            res <- format(res, diff = 1e5)
            if(length(res))
                message(paste(format(res), collapse = "\n"))
        }
    }
    do_exit()
}

.Stangle <- function(args = NULL)
{
    options(warn = 1)
    if (is.null(args)) {
        args <- commandArgs(TRUE)
        args <- paste(args, collapse=" ")
        args <- strsplit(args,'nextArg', fixed = TRUE)[[1L]][-1L]
    }

    Usage <- function() {
        cat("Usage: R CMD Stangle file",
            "",
            "A front-end for Stangle",
            "",
            "Options:",
            "  -h, --help     print this help message and exit",
            "  -v, --version  print version info and exit",
            "  --encoding=enc  assume encoding 'enc' for file",
            "  --options=      comma-separated list of Stangle options",
            "",
            "Report bugs at bugs@r-project.org .",
            sep = "\n")
    }
    do_exit <- function(status = 0L)
        q("no", status = status, runLast = FALSE)

    if (!length(args)) {
        Usage()
        do_exit(1L)
    }

    file <- character()
    encoding <- options <- ""
    while(length(args)) {
        a <- args[1L]
        if (a %in% c("-h", "--help")) {
            Usage()
            do_exit()
        }
        else if (a %in% c("-v", "--version")) {
            cat("Stangle front-end: ",
                R.version[["major"]], ".",  R.version[["minor"]],
                " (r", R.version[["svn rev"]], ")\n", sep = "")
            cat("",
                "Copyright (C) 2006-2011 The R Core Team.",
                "This is free software; see the GNU General Public License version 2",
                "or later for copying conditions.  There is NO warranty.",
                sep = "\n")
            do_exit()
        } else if (substr(a, 1, 11) == "--encoding=") {
            encoding <- substr(a, 12, 1000)
        } else if (substr(a, 1, 10) == "--options=") {
            options <- substr(a, 11, 1000)
        } else if (substr(a, 1, 1) == "-") {
            message(gettextf("Warning: unknown option %s", sQuote(a)),
                    domain = NA)
        } else file <- c(file, a)
        args <- args[-1L]
    }
    if(length(file) != 1L) {
        Usage()
        do_exit(1L)
    }
    args <- list(file)
    args <- c(args, encoding = encoding)
    if(nzchar(options)) {
        opts <- eval(parse(text = paste("list(", options, ")")))
        args <- c(args, opts)
    }
    do.call(Stangle, args)
    do_exit()
}
#   File src/library/utils/R/SweaveDrivers.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

RweaveLatex <- function()
{
    list(setup = RweaveLatexSetup,
         runcode = RweaveLatexRuncode,
         writedoc = RweaveLatexWritedoc,
         finish = RweaveLatexFinish,
         checkopts = RweaveLatexOptions)
}

## We definitely do not want '.' in here, to avoid misidentification
## of file extensions.
.SweaveValidFilenameRegexp <- "^[[:alnum:]#+-_]+$"

RweaveLatexSetup <-
    function(file, syntax, output = NULL, quiet = FALSE, debug = FALSE,
             stylepath, ...)
{
    dots <- list(...)
    if (is.null(output)) {
        prefix.string <- basename(sub(syntax$extension, "", file))
        output <- paste(prefix.string, "tex", sep = ".")
    } else prefix.string <- basename(sub("\\.tex$", "", output))

    if (!quiet) cat("Writing to file ", output, "\n",
                   "Processing code chunks with options ...\n", sep = "")
    encoding <- attr(file, "encoding")
    if (encoding %in% c("ASCII", "bytes")) encoding <- ""
    output <- file(output, open = "w", encoding = encoding)

    if (missing(stylepath)) {
        p <- Sys.getenv("SWEAVE_STYLEPATH_DEFAULT")
        stylepath <-
            if (length(p) >= 1L && nzchar(p[1L])) identical(p, "TRUE") else FALSE
    }
    if (stylepath) {
        styfile <- file.path(R.home("share"), "texmf", "tex", "latex", "Sweave")
        if (.Platform$OS.type == "windows")
            styfile <- chartr("\\", "/", styfile)
        if (length(grep(" ", styfile)))
            warning(gettextf("path to %s contains spaces,\n", sQuote(styfile)),
                    gettext("this may cause problems when running LaTeX"),
                    domain = NA)
    } else styfile <- "Sweave"

    options <- list(prefix = TRUE, prefix.string = prefix.string,
                    engine = "R", print = FALSE, eval = TRUE, fig = FALSE,
                    pdf = TRUE, eps = FALSE, png = FALSE, jpeg = FALSE,
                    grdevice = "", width = 6, height = 6, resolution = 300,
                    term = TRUE, echo = TRUE, keep.source = TRUE,
                    results = "verbatim",
                    split = FALSE, strip.white = "true", include = TRUE,
                    pdf.version = grDevices::pdf.options()$version,
                    pdf.encoding = grDevices::pdf.options()$encoding,
                    pdf.compress = grDevices::pdf.options()$compress,
                    expand = TRUE, # unused by us, for 'highlight'
                    concordance = FALSE, figs.only = TRUE)
    options$.defaults <- options
    options[names(dots)] <- dots

    ## to be on the safe side: see if defaults pass the check
    options <- RweaveLatexOptions(options)

    list(output = output, styfile = styfile, havesty = FALSE,
         haveconcordance = FALSE, debug = debug, quiet = quiet,
         syntax = syntax, options = options,
         chunkout = list(), # a list of open connections
         srclines = integer())
}

makeRweaveLatexCodeRunner <- function(evalFunc = RweaveEvalWithOpt)
{
    ## Return a function suitable as the 'runcode' element
    ## of an Sweave driver.  evalFunc will be used for the
    ## actual evaluation of chunk code.
    ## FIXME: well, actually not for the figures.
    ## If there were just one figure option set, we could eval the chunk
    ## only once.
    function(object, chunk, options) {
        pdf.Swd <- function(name, width, height, ...)
            grDevices::pdf(file = paste(chunkprefix, "pdf", sep = "."),
                           width = width, height = height,
                           version = options$pdf.version,
                           encoding = options$pdf.encoding,
                           compress = options$pdf.compress)
        eps.Swd <- function(name, width, height, ...)
            grDevices::postscript(file = paste(name, "eps", sep = "."),
                                  width = width, height = height,
                                  paper = "special", horizontal = FALSE)
        png.Swd <- function(name, width, height, options, ...)
            grDevices::png(filename = paste(chunkprefix, "png", sep = "."),
                           width = width, height = height,
                           res = options$resolution, units = "in")
        jpeg.Swd <- function(name, width, height, options, ...)
            grDevices::jpeg(filename = paste(chunkprefix, "jpeg", sep = "."),
                            width = width, height = height,
                            res = options$resolution, units = "in")

        if (!(options$engine %in% c("R", "S"))) return(object)

        devs <- devoffs <- list()
        if (options$fig && options$eval) {
            if (options$pdf) {
                devs <- c(devs, list(pdf.Swd))
                devoffs <- c(devoffs, list(grDevices::dev.off))
            }
            if (options$eps) {
                devs <- c(devs, list(eps.Swd))
                devoffs <- c(devoffs, list(grDevices::dev.off))
            }
            if (options$png) {
                devs <- c(devs, list(png.Swd))
                devoffs <- c(devoffs, list(grDevices::dev.off))
            }
            if (options$jpeg) {
                devs <- c(devs, list(jpeg.Swd))
                devoffs <- c(devoffs, list(grDevices::dev.off))
            }
            if (nzchar(grd <- options$grdevice)) {
                devs <- c(devs, list(get(grd, envir = .GlobalEnv)))
                grdo <- paste(grd, "off", sep = ".")
                devoffs <- c(devoffs,
                             if (exists(grdo, envir = .GlobalEnv))
                                 list(get(grdo, envir = .GlobalEnv))
                             else list(grDevices::dev.off))
            }
        }
        if (!object$quiet) {
            cat(formatC(options$chunknr, width = 2), ":")
            if (options$echo) cat(" echo")
            if (options$keep.source) cat(" keep.source")
            if (options$eval) {
                if (options$print) cat(" print")
                if (options$term) cat(" term")
                cat("", options$results)
                if (options$fig) {
                    if (options$eps) cat(" eps")
                    if (options$pdf) cat(" pdf")
                    if (options$png) cat(" png")
                    if (options$jpeg) cat(" jpeg")
                    if (!is.null(options$grdevice)) cat("", options$grdevice)
                }
            }
            cat(" (")
            if (!is.null(options$label))
                cat("label = ", options$label, ", ", sep = "")
            filenum <- attr(chunk, "srcFilenum")[1]
            filename <- attr(chunk, "srcFilenames")[filenum]
            cat(basename(filename), ":", attr(chunk, "srclines")[1], ")", sep = "")
            cat("\n")
        }

        chunkprefix <- RweaveChunkPrefix(options)

        if (options$split) {
            ## [x][[1L]] avoids partial matching of x
            chunkout <- object$chunkout[chunkprefix][[1L]]
            if (is.null(chunkout)) {
                chunkout <- file(paste(chunkprefix, "tex", sep = "."), "w")
                if (!is.null(options$label))
                    object$chunkout[[chunkprefix]] <- chunkout
                if(!grepl(.SweaveValidFilenameRegexp, chunkout))
                    warning("file name ", sQuote(chunkout), " is not portable",
                            call. = FALSE, domain = NA)
            }
        } else chunkout <- object$output

        srcfile <- srcfilecopy(object$filename, chunk, isFile = TRUE)

        ## Note that we edit the error message below, so change both
        ## if you change this line:
        chunkexps <- try(parse(text = chunk, srcfile = srcfile), silent = TRUE)

        if (inherits(chunkexps, "try-error"))
            chunkexps[1L] <- sub(" parse(text = chunk, srcfile = srcfile) : \n ",
                                 "", chunkexps[1L], fixed = TRUE)

        RweaveTryStop(chunkexps, options)

        ## Some worker functions used below...
        putSinput <- function(dce, leading) {
            if (!openSinput) {
                if (!openSchunk) {
                    cat("\\begin{Schunk}\n", file = chunkout)
                    linesout[thisline + 1L] <<- srcline
                    filenumout[thisline + 1L] <<- srcfilenum
                    thisline <<- thisline + 1L
                    openSchunk <<- TRUE
                }
                cat("\\begin{Sinput}", file = chunkout)
                openSinput <<- TRUE
            }
            leading <- max(leading, 1L) # safety check
            cat("\n", paste(getOption("prompt"), dce[seq_len(leading)],
                            sep = "", collapse = "\n"),
                file = chunkout, sep = "")
            if (length(dce) > leading)
                cat("\n", paste(getOption("continue"), dce[-seq_len(leading)],
                                sep = "", collapse = "\n"),
                    file = chunkout, sep = "")
            linesout[thisline + seq_along(dce)] <<- srcline
            filenumout[thisline + seq_along(dce)] <<- srcfilenum
            thisline <<- thisline + length(dce)
        }

        trySrcLines <- function(srcfile, showfrom, showto, ce) {
            lines <- try(suppressWarnings(getSrcLines(srcfile, showfrom, showto)),
                         silent = TRUE)
            if (inherits(lines, "try-error")) {
                if (is.null(ce)) lines <- character()
                else lines <- deparse(ce, width.cutoff = 0.75*getOption("width"))
            }
            lines
        }

        echoComments <- function(showto) {
            if (options$echo && !is.na(lastshown) && lastshown < showto) {
                dce <- trySrcLines(srcfile, lastshown + 1L, showto, NULL)
                linedirs <- grepl("^#line ", dce)
		dce <- dce[!linedirs]
		if (length(dce))
                    putSinput(dce, length(dce)) # These are all trailing comments
                lastshown <<- showto
            }
        }

        openSinput <- FALSE
        openSchunk <- FALSE

        srclines <- attr(chunk, "srclines")
        srcfilenums <- attr(chunk, "srcFilenum")
        linesout <- integer()      # maintains concordance
        filenumout <- integer()	   # ditto
        srcline <- srclines[1L]    # current input line
        srcfilenum <- srcfilenums[1L] # from this file
        thisline <- 0L             # current output line
        lastshown <- 0L            # last line already displayed;

        refline <- NA    # line containing the current named chunk ref
        leading <- 1L    # How many lines get the user prompt

        srcrefs <- attr(chunkexps, "srcref")

        if (length(devs)) {
            if(!grepl(.SweaveValidFilenameRegexp, chunkprefix))
                warning("file name ", sQuote(chunkprefix), " is not portable",
                        call. = FALSE, domain = NA)
            if (options$figs.only)
                devs[[1L]](name = chunkprefix,
                           width = options$width, height = options$height,
                           options)
        }
        SweaveHooks(options, run = TRUE)

        for (nce in seq_along(chunkexps)) {
            ce <- chunkexps[[nce]]
            if (options$keep.source && nce <= length(srcrefs) &&
                !is.null(srcref <- srcrefs[[nce]])) {
                showfrom <- srcref[7L]
                showto <- srcref[8L]

                dce <- trySrcLines(srcfile, lastshown+1L, showto, ce)
                leading <- showfrom - lastshown

                lastshown <- showto
                srcline <- srcref[3L]

                linedirs <- grepl("^#line ", dce)
                dce <- dce[!linedirs]
                # Need to reduce leading lines if some were just removed
                leading <- leading - sum(linedirs[seq_len(leading)])

                while (length(dce) && length(grep("^[[:blank:]]*$", dce[1L]))) {
                    dce <- dce[-1L]
                    leading <- leading - 1L
                }
            } else {
                dce <- deparse(ce, width.cutoff = 0.75*getOption("width"))
                leading <- 1L
            }
            if (object$debug)
                cat("\nRnw> ", paste(dce, collapse = "\n+  "),"\n")

            if (options$echo && length(dce)) putSinput(dce, leading)

            ## avoid the limitations (and overhead) of output text connections
            if (options$eval) {
                tmpcon <- file()
                sink(file = tmpcon)
                err <- evalFunc(ce, options)
                cat("\n")           # make sure final line is complete
                sink()
                output <- readLines(tmpcon)
                close(tmpcon)
                ## delete empty output
                if (length(output) == 1L && !nzchar(output[1L])) output <- NULL
                RweaveTryStop(err, options)
            } else output <- NULL

            ## or writeLines(output)
            if (length(output) && object$debug)
                cat(paste(output, collapse = "\n"))

            if (length(output) && (options$results != "hide")) {
                if (openSinput) {
                    cat("\n\\end{Sinput}\n", file = chunkout)
                    linesout[thisline + 1L:2L] <- srcline
                    filenumout[thisline + 1L:2L] <- srcfilenum
                    thisline <- thisline + 2L
                    openSinput <- FALSE
                }
                if (options$results == "verbatim") {
                    if (!openSchunk) {
                        cat("\\begin{Schunk}\n", file = chunkout)
                        linesout[thisline + 1L] <- srcline
                        filenumout[thisline + 1L] <- srcfilenum
                        thisline <- thisline + 1L
                        openSchunk <- TRUE
                    }
                    cat("\\begin{Soutput}\n", file = chunkout)
                    linesout[thisline + 1L] <- srcline
                    filenumout[thisline + 1L] <- srcfilenum
                    thisline <- thisline + 1L
                }

                output <- paste(output, collapse = "\n")
                if (options$strip.white %in% c("all", "true")) {
                    output <- sub("^[[:space:]]*\n", "", output)
                    output <- sub("\n[[:space:]]*$", "", output)
                    if (options$strip.white == "all")
                        output <- sub("\n[[:space:]]*\n", "\n", output)
                }
                cat(output, file = chunkout)
                count <- sum(strsplit(output, NULL)[[1L]] == "\n")
                if (count > 0L) {
                    linesout[thisline + 1L:count] <- srcline
                    filenumout[thisline + 1L:count] <- srcfilenum
                    thisline <- thisline + count
                }

                remove(output)

                if (options$results == "verbatim") {
                    cat("\n\\end{Soutput}\n", file = chunkout)
                    linesout[thisline + 1L:2L] <- srcline
                    filenumout[thisline + 1L:2L] <- srcfilenum
                    thisline <- thisline + 2L
                }
            }
        } # end of loop over chunkexps.

        ## Echo remaining comments if necessary
        if (options$keep.source) echoComments(length(srcfile$lines))

        if (openSinput) {
            cat("\n\\end{Sinput}\n", file = chunkout)
            linesout[thisline + 1L:2L] <- srcline
            filenumout[thisline + 1L:2L] <- srcfilenum
            thisline <- thisline + 2L
        }

        if (openSchunk) {
            cat("\\end{Schunk}\n", file = chunkout)
            linesout[thisline + 1L] <- srcline
            filenumout[thisline + 1L] <- srcfilenum
            thisline <- thisline + 1L
        }

        if (is.null(options$label) && options$split) close(chunkout)

        if (options$split && options$include) {
            cat("\\input{", chunkprefix, "}\n", sep = "", file = object$output)
            linesout[thisline + 1L] <- srcline
            filenumout[thisline + 1L] <- srcfilenum
            thisline <- thisline + 1L
        }

        if (length(devs)) {
            if (options$figs.only) devoffs[[1L]]()
            for (i in seq_along(devs)) {
                if (options$figs.only && i == 1) next
                devs[[i]](name = chunkprefix, width = options$width,
                          height = options$height, options)
                err <- tryCatch({
                    SweaveHooks(options, run = TRUE)
                    eval(chunkexps, envir = .GlobalEnv)
                }, error = function(e) {
                    devoffs[[i]]()
                    stop(conditionMessage(e), call. = FALSE, domain = NA)
                })
                devoffs[[i]]()
            }

            if (options$include) {
                cat("\\includegraphics{", chunkprefix, "}\n", sep = "",
                    file = object$output)
                linesout[thisline + 1L] <- srcline
                filenumout[thisline + 1L] <- srcfilenum
                thisline <- thisline + 1L
            }
        }
        object$linesout <- c(object$linesout, linesout)
        object$filenumout <- c(object$filenumout, filenumout)
        object
    }
}

RweaveLatexRuncode <- makeRweaveLatexCodeRunner()

RweaveLatexWritedoc <- function(object, chunk)
{
    linesout <- attr(chunk, "srclines")
    filenumout <- attr(chunk, "srcFilenum")

    if (length(grep("\\usepackage[^\\}]*Sweave.*\\}", chunk)))
        object$havesty <- TRUE

    if (!object$havesty) {
 	begindoc <- "^[[:space:]]*\\\\begin\\{document\\}"
 	which <- grep(begindoc, chunk)
 	if (length(which)) {
            chunk[which] <- sub(begindoc,
                                paste("\\\\usepackage{",
                                      object$styfile,
                                      "}\n\\\\begin{document}", sep = ""),
                                chunk[which])
            idx <- c(1L:which, which, seq(from = which+1L,
                     length.out = length(linesout)-which))
            linesout <- linesout[idx]
            filenumout <- filenumout[idx]
            object$havesty <- TRUE
        }
    }

    while(length(pos <- grep(object$syntax$docexpr, chunk)))
    {
        cmdloc <- regexpr(object$syntax$docexpr, chunk[pos[1L]])
        cmd <- substr(chunk[pos[1L]], cmdloc,
                      cmdloc + attr(cmdloc, "match.length") - 1L)
        cmd <- sub(object$syntax$docexpr, "\\1", cmd)
        if (object$options$eval) {
            val <- as.character(eval(parse(text = cmd), envir = .GlobalEnv))
            ## protect against character(), because sub() will fail
            if (length(val) == 0L) val <- ""
        }
        else val <- paste0("\\\\verb#<<", cmd, ">>#")
        ## it's always debatable what \verb delim-character to use;
        ## originally had '{' but that really can mess up LaTeX

        chunk[pos[1L]] <- sub(object$syntax$docexpr, val, chunk[pos[1L]])
    }

    ## Process \SweaveOpts{} or similar
    ## Since they are only supposed to affect code chunks, it is OK
    ## to process all such in a doc chunk at once.
    while(length(pos <- grep(object$syntax$docopt, chunk)))
    {
        opts <- sub(paste0(".*", object$syntax$docopt, ".*"),
                    "\\1", chunk[pos[1L]])
        object$options <- SweaveParseOptions(opts, object$options,
                                             RweaveLatexOptions)

        if (isTRUE(object$options$concordance)
            && !object$haveconcordance) {
            savelabel <- object$options$label
            object$options$label <- "concordance"
            prefix <- RweaveChunkPrefix(object$options)
            object$options$label <- savelabel
            object$concordfile <- paste(prefix, "tex", sep = ".")
            chunk[pos[1L]] <- sub(object$syntax$docopt,
                                  paste0("\\\\input{", prefix, "}"),
                                  chunk[pos[1L]])
            object$haveconcordance <- TRUE
        } else
            chunk[pos[1L]] <- sub(object$syntax$docopt, "", chunk[pos[1L]])
    }

    cat(chunk, sep = "\n", file = object$output)
    object$linesout <- c(object$linesout, linesout)
    object$filenumout <- c(object$filenumout, filenumout)

    object
}

RweaveLatexFinish <- function(object, error = FALSE)
{
    outputname <- summary(object$output)$description
    if (!object$quiet && !error) {
	if(!file.exists(outputname))
	    stop(gettextf("the output file '%s' has disappeared", outputname))
	cat("\n",
	    sprintf("You can now run (pdf)latex on %s", sQuote(outputname)),
	    "\n", sep = "")
    }
    close(object$output)
    if (length(object$chunkout))
        for (con in object$chunkout) close(con)
    if (object$haveconcordance) {
    	## This output format is subject to change.  Currently it contains
    	## three or four parts, separated by colons:
    	## 1.  The output .tex filename
    	## 2.  The input .Rnw filename
    	## 3.  Optionally, the starting line number of the output coded as "ofs nn",
    	##     where nn is the offset to the first output line.  This is omitted if nn is 0.
    	## 4.  The input line numbers corresponding to each output line.
    	##     This are compressed using the following simple scheme:
    	##     The first line number, followed by
    	##     a run-length encoded diff of the rest of the line numbers.
        linesout <- object$linesout
        filenumout <- object$filenumout
        filenames <- object$srcFilenames[filenumout]
        filegps <- rle(filenames)
        offset <- 0L
        for (i in seq_along(filegps$lengths)) {
            len <- filegps$lengths[i]
            inputname <- filegps$values[i]
            vals <- rle(diff(linesout[offset + seq_len(len)]))
            vals <- c(linesout[offset + 1L], as.numeric(rbind(vals$lengths, vals$values)))
    	    concordance <- paste(strwrap(paste(vals, collapse = " ")), collapse = " %\n")
    	    special <- paste0("\\Sconcordance{concordance:", outputname, ":",
                         inputname, ":",
                         if (offset) paste0("ofs ", offset, ":") else "",
                         "%\n",
                         concordance,"}\n")
    	    cat(special, file = object$concordfile, append=offset > 0L)
    	    offset <- offset + len
    	}
    }
    invisible(outputname)
}

## This is the check function for both RweaveLatex and Rtangle drivers
RweaveLatexOptions <- function(options)
{
    defaults <- options[[".defaults"]]

    ## convert a character string to logical
    c2l <- function(x)
        if (is.null(x)) FALSE else suppressWarnings(as.logical(x))

    ## numeric
    NUMOPTS <- c("width", "height", "resolution")

    ## character: largely for safety, but 'label' matters as there
    ## is no default (and someone uses "F")
    CHAROPTS <- c("results", "prefix.string", "engine", "label",
                  "strip.white", "pdf.version", "pdf.encoding", "grdevice")


    for (opt in names(options)) {
        if(opt == ".defaults") next
        oldval <- options[[opt]]
        defval <- defaults[[opt]]
        if(opt %in% CHAROPTS || is.character(defval)) {
        } else if(is.logical(defval))
            options[[opt]] <- c2l(oldval)
        else if(opt %in% NUMOPTS || is.numeric(defval))
            options[[opt]] <- as.numeric(oldval)
        else if(!is.na(newval <- c2l(oldval)))
            options[[opt]] <- newval
        else if(!is.na(newval <- suppressWarnings(as.numeric(oldval))))
            options[[opt]] <- newval
        if (is.na(options[[opt]]))
            stop(gettextf("invalid value for %s : %s", sQuote(opt), oldval),
                 domain = NA)
    }

    if (!is.null(options$results)) {
        res <- as.character(options$results)
        if(tolower(res) != res) # documented as lower-case
            warning("value of 'results' option should be lowercase",
                    call. = FALSE)
        options$results <- tolower(res)
    }
    options$results <- match.arg(options$results, c("verbatim", "tex", "hide"))

    if (!is.null(options$strip.white)) {
        res <- as.character(options$strip.white)
        if(tolower(res) != res)
            warning("value of 'strip.white' option should be lowercase",
                    call. = FALSE)
        options$strip.white <- tolower(res)
    }
    options$strip.white <-
        match.arg(options$strip.white, c("true", "false", "all"))
    options
}


RweaveChunkPrefix <- function(options)
{
    if (!is.null(options$label)) {
        if (options$prefix)
            chunkprefix <- paste0(options$prefix.string, "-", options$label)
        else
            chunkprefix <- options$label
    } else
        chunkprefix <- paste0(options$prefix.string, "-",
                              formatC(options$chunknr, flag = "0", width = 3))
    chunkprefix
}

RweaveEvalWithOpt <- function (expr, options)
{
    if (options$eval) {
        res <- try(withVisible(eval(expr, .GlobalEnv)), silent = TRUE)
        if (inherits(res, "try-error")) return(res)
        if (options$print || (options$term && res$visible)) {
            if (.isMethodsDispatchOn() && isS4(res$value))
                methods:::show(res$value) else print(res$value)
        }
    }
    res
}

RweaveTryStop <- function(err, options)
{
    if (inherits(err, "try-error")) {
        cat("\n")
        msg <- paste(" chunk", options$chunknr)
        if (!is.null(options$label))
            msg <- paste0(msg, " (label = ", options$label, ")")
        msg <- paste(msg, "\n")
        stop(msg, err, call. = FALSE)
    }
}

###------------------------------------------------------------------------

Rtangle <-  function()
{
    list(setup = RtangleSetup,
         runcode = RtangleRuncode,
         writedoc = RtangleWritedoc,
         finish = RtangleFinish,
         checkopts = RweaveLatexOptions)
}


RtangleSetup <-
    function(file, syntax, output = NULL, annotate = TRUE, split = FALSE,
             quiet = FALSE, ...)
{
    dots <- list(...)
    if (is.null(output)) {
        prefix.string <- basename(sub(syntax$extension, "", file))
        ## This is odd, since for split = TRUE it uses the engine name.
        output <- paste(prefix.string, "R", sep = ".")
    } else
        prefix.string <- basename(sub("\\.[rsRS]$", "", output))

    if (!split) {
        if (identical(output, "stdout")) output <- stdout()
        else if (identical(output, "stderr")) output <- stderr()
        else {
            if (!quiet) cat("Writing to file", output, "\n")
            ## We could at some future point try to write the file in
            ## 'encoding'.
            output <- file(output, open = "w")
        }
        lines <- c(sprintf("R code from vignette source '%s'", file),
                   if(attr(file, "encoding") != "ASCII")
                   sprintf("Encoding: %s", localeToCharset()[1L])
                   )
        lines <- c(paste("###", lines), "")
        writeLines(lines, output)
    } else {
        if (!quiet) cat("Writing chunks to files ...\n")
        output <- NULL
    }

    options <- list(split = split, prefix = TRUE,
                    prefix.string = prefix.string,
                    engine = "R", eval = TRUE,
                    show.line.nos = FALSE)
    options$.defaults <- options
    options[names(dots)] <- dots

    ## to be on the safe side: see if defaults pass the check
    options <- RweaveLatexOptions(options)

    list(output = output, annotate = annotate, options = options,
         chunkout = list(), quiet = quiet, syntax = syntax)
}


RtangleRuncode <-  function(object, chunk, options)
{
    if (!(options$engine %in% c("R", "S"))) return(object)

    chunkprefix <- RweaveChunkPrefix(options)

    if (options$split) {
        if(!grepl(.SweaveValidFilenameRegexp, chunkprefix))
            warning("file name ", sQuote(chunkprefix), " is not portable",
                    call. = FALSE, domain = NA)
        outfile <- paste(chunkprefix, options$engine, sep = ".")
        if (!object$quiet) cat(options$chunknr, ":", outfile,"\n")
        ## [x][[1L]] avoids partial matching of x
        chunkout <- object$chunkout[chunkprefix][[1L]]
        if (is.null(chunkout)) {
            chunkout <- file(outfile, "w")
            if (!is.null(options$label))
                object$chunkout[[chunkprefix]] <- chunkout
        }
    } else
        chunkout <- object$output

    if (object$annotate) {
        lnos <- grep("^#line ", chunk, value = TRUE)
        if(length(lnos)) {
            srclines <- attr(chunk, "srclines")
            srcfilenum <- attr(chunk, "srcFilenum")
            ## this currently includes the chunk header
            lno <- if (length(srclines)) paste(min(srclines), max(srclines), sep = "-") else srclines
            fn <- sub('[^"]*"([^"]+).*', "\\1", lnos[1L])
        }
        cat("###################################################\n",
            "### code chunk number ", options$chunknr,
            ": ",
            if(!is.null(options$label)) options$label
            else paste(fn, lno, sep = ":"),
            ifelse(options$eval, "", " (eval = FALSE)"), "\n",
            "###################################################\n",
            file = chunkout, sep = "")
    }

    ## The next returns a character vector of the logical options
    ## which are true and have hooks set.
    hooks <- SweaveHooks(options, run = FALSE)
    for (k in hooks)
        cat("getOption(\"SweaveHooks\")[[\"", k, "\"]]()\n",
            file = chunkout, sep = "")

    if (!options$show.line.nos)
        chunk <- grep("^#line ", chunk, value = TRUE, invert = TRUE)
    if (!options$eval) chunk <- paste("##", chunk)
    cat(chunk, "\n", file = chunkout, sep = "\n")
    if (is.null(options$label) && options$split) close(chunkout)
    object
}

RtangleWritedoc <- function(object, chunk)
{
    while(length(pos <- grep(object$syntax$docopt, chunk))) {
        opts <- sub(paste0(".*", object$syntax$docopt, ".*"),
                    "\\1", chunk[pos[1L]])
        object$options <- SweaveParseOptions(opts, object$options,
                                             RweaveLatexOptions)
        chunk[pos[1L]] <- sub(object$syntax$docopt, "", chunk[pos[1L]])
    }
    object
}


RtangleFinish <- function(object, error = FALSE)
{
    ## might be stdout() or stderr()
    if (!is.null(object$output) && object$output >= 3)
        close(object$output)

    if (length(object$chunkout))
        for (con in object$chunkout) close(con)
}
#  File src/library/utils/R/URLencode.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

URLencode <- function(URL, reserved = FALSE)
{
    ## It is unsafe to use ranges here as collation is locale-dependent.
    ## We want to do this on characters and not on bytes.
    OK <- paste0("[^-ABCDEFGHIJKLMNOPQRSTUVWXYZ",
		"abcdefghijklmnopqrstuvwxyz0123456789$_.+!*'(),",
		if(!reserved) ";/?:@=&", "]")
    x <- strsplit(URL, "")[[1L]]
    z <- grep(OK, x)
    if(length(z)) {
        y <- sapply(x[z], function(x)
                    paste0("%", as.character(charToRaw(x)), collapse = ""))
        x[z] <- y
    }
    paste(x, collapse="")
}

URLdecode <- function(URL)
{
    x <- charToRaw(URL)
    pc <- charToRaw("%")
    out <- raw(0L)
    i <- 1L
    while(i <= length(x)) {
        if(x[i] != pc) {
            out <- c(out, x[i])
            i <- i + 1L
        } else {
            y <- as.integer(x[i + 1L:2L])
            y[y > 96L] <- y[y > 96L] - 32L # a-f -> A-F
            y[y > 57L] <- y[y > 57L] - 7L  # A-F
            y <- sum((y - 48L) * c(16L, 1L))
            out <- c(out, as.raw(as.character(y)))
            i <- i + 3L
        }
    }
    rawToChar(out)
}
#  File src/library/utils/R/adist.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

adist <-
function(x, y = NULL, costs = NULL, counts = FALSE, fixed = TRUE,
         partial = !fixed, ignore.case = FALSE, useBytes = FALSE)
{
    bytesToInt <- function(x) {
        if(is.na(x)) return(NA_integer_)
        as.integer(charToRaw(x))
    }

    costs <- .amatch_costs(costs)

    nmx <- names(x)
    x <- as.character(x)
    names(x) <- nmx

    if(!is.null(y)) {
        nmy <- names(y)
        y <- as.character(y)
        names(y) <- nmy
    }

    if(!identical(fixed, FALSE) && !identical(partial, TRUE)) {
        ex <- Encoding(x)
        useBytes <- identical(useBytes, TRUE) || any(ex == "bytes")
        if(!is.null(y)) {
            ey <- Encoding(y)
            useBytes <- useBytes || any(ey == "bytes")
        }
        if(useBytes) {
            x <- lapply(x, bytesToInt)
            y <- if(is.null(y)) {
                x
            } else {
                lapply(y, bytesToInt)
            }
        } else {
            ignore.case <- identical(ignore.case, TRUE)
            x <- if(ignore.case) {
                lapply(tolower(enc2utf8(x)), utf8ToInt)
            } else {
                lapply(enc2utf8(x), utf8ToInt)
            }
            y <- if(is.null(y)) {
                x
            } else if(ignore.case) {
                lapply(tolower(enc2utf8(y)), utf8ToInt)
            } else {
                lapply(enc2utf8(y), utf8ToInt)
            }
        }
    }
    else {
        if(is.null(y)) {
            y <- x
        }
        ## TRE needs integer costs: coerce here for simplicity.
        costs <- as.integer(costs)
    }

    .Internal(adist(x, y, costs, counts, fixed, partial, ignore.case,
                    useBytes))
}

aregexec <-
function(pattern, text, max.distance = 0.1, costs = NULL,
         ignore.case = FALSE, fixed = FALSE, useBytes = FALSE)
{
    ## TRE needs integer costs: coerce here for simplicity.
    costs <- as.integer(.amatch_costs(costs))
    bounds <- .amatch_bounds(max.distance)

    .Internal(aregexec(as.character(pattern),
                       as.character(text),
                       bounds, costs, ignore.case, fixed, useBytes))
}

## No longer used by adist(), but could be more generally useful ...

regquote <-
function(x)
    gsub("([*.?+^&\\[])", "\\\\\\1", x)

#  File src/library/utils/R/alarm.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

alarm <- function() {cat("\a"); flush.console()}

#  File src/library/utils/R/apropos.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

apropos <- function (what, where = FALSE, ignore.case = TRUE, mode = "any")
{
    stopifnot(is.character(what))
    x <- character(0L)
    check.mode <- mode != "any"
    for (i in seq_along(search())) {
	li <-
	    if(ignore.case)
		grep(what, ls(pos = i, all.names = TRUE),
		     ignore.case = TRUE, value = TRUE)
	    else ls(pos = i, pattern = what, all.names = TRUE)
	if(length(li)) {
	    if(check.mode)
		li <- li[sapply(li, exists, where = i,
				mode = mode, inherits = FALSE)]
	    x <- c(x, if(where) structure(li, names = rep.int(i, length(li))) else li)
	}
    }
    sort(x)
}

find <- function(what, mode = "any", numeric = FALSE, simple.words=TRUE)
{
    stopifnot(is.character(what))
    if(length(what) > 1L) {
        warning("elements of 'what' after the first will be ignored")
        what <- what[1L]
    }
#   would need to escape at least + * | as well
#     if(simple.words)
# 	what <- gsub("([.[])", "\\\\\\1", paste0("^",what,"$"))
    len.s <- length(sp <- search())
    ind <- logical(len.s)
    check.mode <- mode != "any"
    for (i in 1L:len.s) {
        if(simple.words) {
            found <- what %in% ls(pos = i, all.names = TRUE)
            if(found && check.mode)
                found <- exists(what, where = i, mode = mode, inherits=FALSE)
            ind[i] <- found
        } else {
            li <- ls(pos = i, pattern = what, all.names = TRUE)
            ll <- length(li)
            if(ll > 0 && check.mode) {
                mode.ok <- sapply(li, exists, where = i, mode = mode,
                                  inherits = FALSE)
                ll <- sum(mode.ok)
                if(ll >= 2) # some languages have multiple plurals
                    warning(sprintf(ngettext(ll,
                                             "%d occurrence in %s",
                                             "%d occurrences in %s"), ll, sp[i]),
                            domain = NA)
            }
            ind[i] <- ll > 0L
        }
    }
    ## found name in  search()[ ind ]
    if(numeric) structure(which(ind), names=sp[ind]) else sp[ind]
}

#  File src/library/utils/R/aspell.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


aspell <-
function(files, filter, control = list(), encoding = "unknown",
         program = NULL, dictionaries = character())
{
    ## Take the given files and feed them through spell checker in
    ## Ispell pipe mode.

    ## Think about options and more command line options eventually.

    program <- aspell_find_program(program)
    if(is.na(program))
        stop("No suitable spell-checker program found")

    ## Be nice.
    if(inherits(files, "Rd"))
        files <- list(files)

    files_are_names <- is.character(files)

    filter_args <- list()
    if(missing(filter) || is.null(filter)) {
        filter <- if(!files_are_names) {
            function(ifile, encoding) {
                if(inherits(ifile, "srcfile"))
                    readLines(ifile$filename, encoding = encoding,
                              warn = FALSE)
                else if(inherits(ifile, "connection"))
                    readLines(ifile, encoding = encoding)
                else {
                    ## What should this do with encodings?
                    as.character(ifile)
                }
            }
        }
        else NULL
    }
    else if(is.character(filter)) {
        ## Look up filter in aspell filter db.
        filter_name <- filter[1L]
        filter <- aspell_filter_db[[filter_name]]
        ## Warn if the filter was not found in the db.
        if(is.null(filter))
            warning(gettextf("Filter '%s' is not available.",
                             filter_name),
                    domain = NA)
    }
    else if(is.list(filter)) {
        ## Support
        ##   list("Rd", drop = "\\references"
        ## at least for now.
        filter_name <- filter[[1L]][1L]
        filter_args <- filter[-1L]
        filter <- aspell_filter_db[[filter_name]]
        ## Warn if the filter was not found in the db.
        if(is.null(filter))
            warning(gettextf("Filter '%s' is not available.",
                             filter_name),
                    domain = NA)
    }
    else if(!is.function(filter))
        stop("Invalid 'filter' argument.")

    encoding <- rep(encoding, length.out = length(files))

    verbose <- getOption("verbose")

    db <- data.frame(Original = character(), File = character(),
                     Line = integer(), Column = integer(),
                     stringsAsFactors = FALSE)
    db$Suggestions <- list()

    tfile <- tempfile("aspell")
    on.exit(unlink(tfile))

    if(length(dictionaries)) {
        paths <- aspell_find_dictionaries(dictionaries)
        ind <- paths == ""
        if(any(ind)) {
            warning(gettextf("The following dictionaries were not found:\n%s",
                             paste(sprintf("  %s", dictionaries[ind]),
                                   collapse = "\n")),
                    domain = NA)
            paths <- paths[!ind]
        }
        if(length(paths)) {
            words <- unlist(lapply(paths, readRDS), use.names = FALSE)
            personal <- tempfile("aspell_personal")
            on.exit(unlink(personal), add = TRUE)
            ## <FIXME>
            ## How can we get the right language set (if needed)?
            ## Maybe aspell() needs an additional 'language' arg?
            aspell_write_personal_dictionary_file(words, personal,
                                                  program = program)
            ## </FIXME>
            control <- c(control, "-p", shQuote(personal))
        }
    }

    ## No special expansion of control argument for now.
    control <- as.character(control)

    fnames <- names(files)
    files <- as.list(files)

    for (i in seq_along(files)) {

        file <- files[[i]]
        if(files_are_names)
            fname <- file
        else {
            ## Try srcfiles and srcrefs ...
            fname <- if(inherits(file, "srcfile"))
                file$filename
            else
                attr(attr(file, "srcref"), "srcfile")$filename
            ## As a last resort, try the names of the files argument.
            if(is.null(fname))
                fname <- fnames[i]
            ## If unknown ...
            if(is.null(fname))
                fname <- "<unknown>"
        }

        enc <- encoding[i]

        if(verbose)
            message(gettextf("Processing file %s", fname),
                    domain = NA)

        lines <- if(is.null(filter))
            readLines(file, encoding = enc)
        else {
            ## Assume that filter takes an input file (and additional
            ## arguments) and return a character vector.
            do.call(filter, c(list(file, encoding = enc), filter_args))
        }

        ## Need to escape all lines with carets to ensure Aspell handles
        ## them as data: the Aspell docs say
        ##   It is recommended that programmatic interfaces prefix every
        ##   data line with an uparrow to protect themselves against
        ##   future changes in Aspell.
        writeLines(paste0("^", lines), tfile)
        ## Note that this re-encodes character strings with marked
        ## encodings to the current encoding (which is definitely fine
        ## if this is UTF-8 and Aspell was compiled with full UTF-8
        ## support).  Alternatively, we could try using something along
        ## the lines of
        ##   writeLines(paste0("^", lines), tfile,
        ##              useBytes = TRUE)
        ## and pass the encoding info to Aspell in case we know it.

        out <- tools:::.system_with_capture(program, c("-a", control),
                                            stdin = tfile)
                                           
	if(out$status != 0L)
	    stop(gettextf("Running aspell failed with diagnostics:\n%s",
			  paste(out$stderr, collapse = "\n")),
                 domain = NA)

	## Hopefully everything worked ok.
	lines <- out$stdout[-1L]
	pos <- cumsum(lines == "") + 1L

	## Format is as follows.
	## First line is a header.
	## Blank lines separate the results for each line.
	## Results for the word on each line are given as follows.
	## * If the word was found in the main dictionary, or your personal
	##   dictionary, then the line contains only a `*'.
	## * If the word is not in the dictionary, but there are
	##   suggestions, then the line contains an `&', a space, the
	##   misspelled word, a space, the number of near misses, the number
	##   of characters between the beginning of the line and the
	##   beginning of the misspelled word, a colon, another space, and a
	##   list of the suggestions separated by commas and spaces.
	## * If the word does not appear in the dictionary, and there are no
	##   suggestions, then the line contains a `#', a space, the
	##   misspelled word, a space, and the character offset from the
	##   beginning of the line.
	## This can be summarized as follows:
	##   OK: *
	##   Suggestions: & original count offset: miss, miss, ...
	##   None: # original offset

	## Look at words not in dictionary with suggestions.

	ind <- grepl("^&", lines)
	if(any(ind)) {
	    info <- strsplit(lines[ind], ": ", fixed = TRUE)
	    one <- strsplit(sapply(info, `[`, 1L), " ",  fixed = TRUE)
	    two <- strsplit(sapply(info, `[`, 2L), ", ", fixed = TRUE)
	    db1 <- data.frame(Original =
			      as.character(sapply(one, `[`, 2L)),
			      File = fname,
			      Line = pos[ind],
			      Column =
			      as.integer(sapply(one, `[`, 4L)),
			      stringsAsFactors = FALSE)
	    db1$Suggestions <- two
	    db <- rbind(db, db1)
	}
	## Looks at words not in dictionary with no suggestions.
	ind <- grepl("^#", lines)
	if(any(ind)) {
	    one <- strsplit(lines[ind], " ", fixed = TRUE)
	    db1 <- data.frame(Original =
			      as.character(sapply(one, `[`, 2L)),
			      File = fname,
			      Line = pos[ind],
			      Column =
			      as.integer(sapply(one, `[`, 3L)),
			      stringsAsFactors = FALSE)
	    db1$Suggestions <- vector("list", length(one))
	    db <- rbind(db, db1)
	}
    }

    class(db) <- c("aspell", "data.frame")
    db
}

print.aspell <-
function(x, sort = TRUE, verbose = FALSE, indent = 2L, ...)
{
    ## A very simple printer ...
    if(!(nr <- nrow(x))) return(invisible(x))

    if (sort)
    	x <- x[order(x$Original, x$File, x$Line, x$Column), ]

    if (verbose)
    	out <-
    	    sprintf("%sWord: %s (%s:%d:%d)\n%s",
    	            c("", rep.int("\n", nr - 1L)),
    	            x$Original, x$File, x$Line, x$Column,
    	            formatDL(rep.int("Suggestions", nr),
    	                     sapply(x$Suggestions, paste, collapse = " "),
    	                     style = "list"))
    else {
        s <- split(sprintf("%s:%d:%d", x$File, x$Line, x$Column),
                   x$Original)
        sep <- sprintf("\n%s",
                       paste(rep.int(" ", indent), collapse = ""))
        out <- paste(names(s),
                     sapply(s, paste, collapse = sep),
                     sep = sep, collapse = "\n\n")
    }
    writeLines(out)
    invisible(x)
}

summary.aspell <-
function(object, ...)
{
    words <- sort(unique(object$Original))
    if(length(words)) {
        writeLines("Possibly mis-spelled words:")
        print(words)
    }
    invisible(words)
}

aspell_filter_db <- new.env(hash = FALSE) # small
aspell_filter_db$Rd <- tools::RdTextFilter
aspell_filter_db$Sweave <- tools::SweaveTeXFilter

aspell_find_program <-
function(program = NULL)
{
    check <- !is.null(program) || !is.null(names(program))
    if(is.null(program))
        program <- getOption("aspell_program")
    if(is.null(program))
        program <- c("aspell", "hunspell", "ispell")
    program <- Filter(nzchar, Sys.which(program))[1L]
    if(!is.na(program) && check) {
        out <- c(system(sprintf("%s -v", program),
                        intern = TRUE), "")[1L]
        if(grepl("really Aspell", out))
            names(program) <- "aspell"
        else if(grepl("really Hunspell", out))
            names(program) <- "hunspell"
        else if(grepl("International Ispell", out))
            names(program) <- "ispell"
        else
            names(program) <- NA_character_
    }
    program
}

aspell_dictionaries_R <- "en_stats"

aspell_find_dictionaries <-
function(dictionaries, dirnames = character())
{
    dictionaries <- as.character(dictionaries)
    if(!(n <- length(dictionaries))) return(character())

    ## Always search the R system dictionary directory first.
    dirnames <- c(file.path(R.home("share"), "dictionaries"), dirnames)

    ## For now, all dictionary files should be .rds files.
    ind <- !grepl("\\.rds$", dictionaries)
    if(any(ind))
        dictionaries[ind] <- sprintf("%s.rds", dictionaries[ind])

    out <- character(n)
    ## Dictionaries with no path separators are looked for in the given
    ## dictionary directories (by default, the R system dictionary
    ## directory).
    ind <- grepl(.Platform$file.sep, dictionaries, fixed = TRUE)
    ## (Equivalently, could check where paths == basename(paths).)
    if(length(pos <- which(ind))) {
        pos <- pos[file_test("-f", dictionaries[pos])]
        out[pos] <- normalizePath(dictionaries[pos], "/")
    }
    if(length(pos <- which(!ind))) {
        out[pos] <- find_files_in_directories(dictionaries[pos],
                                              dirnames)
    }

    out
}

### Utilities.

aspell_inspect_context <-
function(x)
{
    x <- split(x, x$File)
    y <- Map(function(f, x) {
        lines <- readLines(f, warn = FALSE)[x$Line]
        cbind(f,
              x$Line,
              substring(lines, 1L, x$Column - 1L),
              x$Original,
              substring(lines, x$Column + nchar(x$Original)))
    },
             names(x), x)
    y <- data.frame(do.call(rbind, y), stringsAsFactors = FALSE)
    names(y) <- c("File", "Line", "Left", "Original", "Right")
    class(y) <- c("aspell_inspect_context", "data.frame")
    y
}

print.aspell_inspect_context <-
function(x, ...)
{
    s <- split(x, x$File)
    nms <- names(s)
    for(i in seq_along(s)) {
        e <- s[[i]]
        writeLines(c(sprintf("File '%s':", nms[i]),
                     sprintf("  Line %s: \"%s\", \"%s\", \"%s\"",
                             format(e$Line),
                             gsub("\"", "\\\"", e$Left),
                             e$Original,
                             gsub("\"", "\\\"", e$Right)),
                     ""))
    }
    invisible(x)
}


## For spell-checking the R manuals:

## This can really only be done with Aspell as the other checkers have
## no texinfo mode.

aspell_control_R_manuals <-
    list(aspell =
         c("--master=en_US",
           "--add-extra-dicts=en_GB",
           "--mode=texinfo",
           "--add-texinfo-ignore=acronym",
           "--add-texinfo-ignore=deftypefun",
           "--add-texinfo-ignore=deftypefunx",
           "--add-texinfo-ignore=findex",
           "--add-texinfo-ignore=enindex",
           "--add-texinfo-ignore=include",
           "--add-texinfo-ignore=ifclear",
           "--add-texinfo-ignore=ifset",
           "--add-texinfo-ignore=math",
           "--add-texinfo-ignore=macro",
           "--add-texinfo-ignore=multitable",
           "--add-texinfo-ignore=node",
           "--add-texinfo-ignore=pkg",
           "--add-texinfo-ignore=printindex",
           "--add-texinfo-ignore=set",
           "--add-texinfo-ignore=vindex",
           "--add-texinfo-ignore-env=menu",
           "--add-texinfo-ignore=CRANpkg"
           ),
         hunspell =
         c("-d en_US,en_GB"))

aspell_R_manuals <-
function(which = NULL, dir = NULL, program = NULL,
         dictionaries = aspell_dictionaries_R)
{
    if(is.null(dir)) dir <- tools:::.R_top_srcdir_from_Rd()
    ## Allow specifying 'R-exts' and alikes, or full paths.
    files <- if(is.null(which)) {
        Sys.glob(file.path(dir, "doc", "manual", "*.texi"))
    } else {
        ind <- which(which ==
                     basename(tools::file_path_sans_ext(which)))
        which[ind] <-
            file.path(dir, "doc", "manual",
                      sprintf("%s.texi", which[ind]))
        which
    }

    program <- aspell_find_program(program)

    aspell(files,
           control = aspell_control_R_manuals[[names(program)]],
           program = program,
           dictionaries = dictionaries)
}

## For spell-checking the R Rd files:

aspell_control_R_Rd_files <-
    list(aspell =
         c("--master=en_US",
           "--add-extra-dicts=en_GB"),
         hunspell =
         c("-d en_US,en_GB"))

aspell_R_Rd_files <-
function(which = NULL, dir = NULL, drop = "\\references",
         program = NULL, dictionaries = aspell_dictionaries_R)
{
    files <- character()

    if(is.null(dir)) dir <- tools:::.R_top_srcdir_from_Rd()

    if(is.null(which)) {
        which <- tools:::.get_standard_package_names()$base
        # CHANGES.Rd could be dropped from checks in the future;
        # it will not be updated post 2.15.0
        files <- c(file.path(dir, "doc", "NEWS.Rd"),
                   file.path(dir, "src", "gnuwin32", "CHANGES.Rd"))
        files <- files[file_test("-f", files)]
    }

    files <-
        c(files,
          unlist(lapply(file.path(dir, "src", "library", which, "man"),
                        tools::list_files_with_type,
                        "docs", OS_subdirs = c("unix", "windows")),
                 use.names = FALSE))

    program <- aspell_find_program(program)

    aspell(files,
           filter = list("Rd", drop = drop),
           control = aspell_control_R_Rd_files[[names(program)]],
           program = program,
           dictionaries = dictionaries)
}

## For spell-checking Rd files in a package:

aspell_package_Rd_files <-
function(dir, drop = c("\\author", "\\references"),
         control = list(), program = NULL, dictionaries = character())
{
    dir <- normalizePath(dir, "/")

    subdir <- file.path(dir, "man")
    files <- if(file_test("-d", subdir))
        tools::list_files_with_type(subdir,
                                    "docs",
                                    OS_subdirs = c("unix", "windows"))
    else character()

    meta <- tools:::.get_package_metadata(dir, installed = FALSE)
    if(is.na(encoding <- meta["Encoding"]))
        encoding <- "unknown"

    defaults <- .aspell_package_defaults(dir, encoding)$Rd_files
    if(!is.null(defaults)) {
        ## Direct settings currently override (could add a list add =
        ## TRUE mechanism eventually).
        if(!is.null(d <- defaults$drop))
            drop <- d
        if(!is.null(d <- defaults$control))
            control <- d
        if(!is.null(d <- defaults$program))
            program <- d
        if(!is.null(d <- defaults$dictionaries)) {
            dictionaries <-
                aspell_find_dictionaries(d, file.path(dir, ".aspell"))
        }
        ## <FIXME>
        ## Deprecated in favor of specifying R level dictionaries.
        ## Maybe give a warning (in particular if both are given)?
        if(!is.null(d <- defaults$personal))
            control <- c(control,
                         sprintf("-p %s",
                                 shQuote(file.path(dir, ".aspell", d))))
        ## </FIXME>
    }

    aspell(files,
           filter = list("Rd", drop = drop),
           control = control,
           encoding = encoding,
           program = program,
           dictionaries = dictionaries)
}

## For spell-checking the R vignettes:

## This should really be done with Aspell as the other checkers have far
## less powerful TeX modes.

aspell_control_R_vignettes <-
    list(aspell =
         c("--mode=tex",
           "--master=en_US",
           "--add-extra-dicts=en_GB",
           "--add-tex-command='code p'",
           "--add-tex-command='pkg p'",
           "--add-tex-command='CRANpkg p'"
           ),
         hunspell =
         c("-t", "-d en_US,en_GB"))

aspell_R_vignettes <-
function(program = NULL, dictionaries = aspell_dictionaries_R)
{
    files <- Sys.glob(file.path(tools:::.R_top_srcdir_from_Rd(),
                                "src", "library", "*", "vignettes",
                                "*.Rnw"))

    program <- aspell_find_program(program)

    aspell(files,
           filter = "Sweave",
           control = aspell_control_R_vignettes[[names(program)]],
           program = program,
           dictionaries = dictionaries)
}

## For spell-checking vignettes in a package:

## This should really be done with Aspell as the other checkers have far
## less powerful TeX modes.

aspell_control_package_vignettes <-
    list(aspell =
         c("--add-tex-command='citep oop'",
           "--add-tex-command='Sexpr p'",
           "--add-tex-command='code p'",
           "--add-tex-command='pkg p'",
           "--add-tex-command='proglang p'",
           "--add-tex-command='samp p'"
           ))

aspell_package_vignettes <-
function(dir,
         control = list(), program = NULL, dictionaries = character())
{
    dir <- tools::file_path_as_absolute(dir)

    subdir <- file.path(dir, "inst", "doc")
    files <- if(file_test("-d", subdir))
        tools::list_files_with_type(subdir, "vignette")
    else character()

    meta <- tools:::.get_package_metadata(dir, installed = FALSE)
    if(is.na(encoding <- meta["Encoding"]))
        encoding <- "unknown"

    defaults <- .aspell_package_defaults(dir, encoding)$vignettes
    if(!is.null(defaults)) {
        if(!is.null(d <- defaults$control))
            control <- d
        if(!is.null(d <- defaults$program))
            program <- d
        if(!is.null(d <- defaults$dictionaries)) {
            dictionaries <-
                aspell_find_dictionaries(d, file.path(dir, ".aspell"))
        }
        ## <FIXME>
        ## Deprecated in favor of specifying R level dictionaries.
        ## Maybe give a warning (in particular if both are given)?
        if(!is.null(d <- defaults$personal))
            control <- c(control,
                         sprintf("-p %s",
                                 shQuote(file.path(dir, ".aspell", d))))
        ## </FIXME>
    }

    program <- aspell_find_program(program)

    aspell(files,
           filter = "Sweave",
           control =
           c("-t",
             aspell_control_package_vignettes[[names(program)]],
             control),
           program = program,
           dictionaries = dictionaries)
}

## Spell-checking R files.

aspell_filter_db$R <-
function(ifile, encoding = "unknown", ignore = character())
{
    pd <- get_parse_data_for_message_strings(ifile, encoding)
    if(is.null(pd) || !NROW(pd)) return(character())

    ## Strip the string delimiters.
    pd$text <- substring(pd$text, 2L, nchar(pd$text) - 1L)
    ## Replace whitespace C backslash escape sequences by whitespace.
    pd$text <- gsub("(^|[^\\])\\\\[fnrt]", "\\1  ", pd$text)
    pd$text <- gsub(  "([^\\])\\\\[fnrt]", "\\1  ", pd$text)
    ## (Do this twice for now because in e.g.
    ##    \n\t\tInformation on package %s
    ## the first \t is not matched the first time.  Alternatively, we
    ## could match with
    ##    (^|[^\\])((\\\\[fnrt])+)
    ## but then computing the replacement (\\1 plus as many blanks as
    ## the characters in \\2) is not straightforward.
    ## For gettextf() calls, replace basic percent escape sequences by
    ## whitespace.
    ind <- pd$caller == "gettextf"
    if(any(ind)) {
        pd$text[ind] <-
            gsub("(^|[^%])%[dioxXfeEgGaAs]", "\\1  ", pd$text[ind])
        pd$text[ind] <-
            gsub("  ([^%])%[dioxXfeEgGaAs]", "\\1  ", pd$text[ind])
        ## (See above for doing this twice.)
    }

    lines <- readLines(ifile, encoding = encoding)

    ## Column positions in the parse data have tabs expanded to tab
    ## stops using a tab width of 8, so for lines with tabs we need to
    ## map the column positions back to character positions.
    lines_in_pd <- sort(unique(c(pd$line1, pd$line2)))
    tab <- Map(function(tp, nc) {
        if(tp[1L] == -1L) return(NULL)
        widths <- rep.int(1, nc)
        for(i in tp) {
            cols <- cumsum(widths)
            widths[i] <- 8 - (cols[i] - 1) %% 8
        }
        cumsum(c(1, widths))
    },
               gregexpr("\t", lines[lines_in_pd], fixed = TRUE),
               nchar(lines[lines_in_pd]))
    names(tab) <- lines_in_pd

    lines[lines_in_pd] <- gsub("[^\t]", " ", lines[lines_in_pd])
    lines[-lines_in_pd] <- ""

    for(entry in split(pd, seq_len(NROW(pd)))) {
        line1 <- entry$line1
        line2 <- entry$line2
        col1 <- entry$col1 + 1L
        col2 <- entry$col2 - 1L
        if(line1 == line2) {
            if(length(ptab <- tab[[as.character(line1)]])) {
                col1 <- which(ptab == col1)
                col2 <- which(ptab == col2)
            }
            substring(lines[line1], col1, col2) <- entry$text
        } else {
            texts <- unlist(strsplit(entry$text, "\n", fixed = TRUE))
            n <- length(texts)
            if(length(ptab <- tab[[as.character(line1)]])) {
                col1 <- which(ptab == col1)
            }
            substring(lines[line1], col1) <- texts[1L]
            pos <- seq(from = 2, length.out = n - 2)
            if(length(pos))
                lines[line1 + pos - 1] <- texts[pos]
            if(length(ptab <- tab[[as.character(line2)]])) {
                col2 <- which(ptab == col2)
            }
            substring(lines[line2], 1L, col2) <- texts[n]
        }
    }

    for(re in ignore[nzchar(ignore)])
        lines <- blank_out_regexp_matches(lines, re)

    lines
}

get_parse_data_for_message_strings <-
function(file, encoding = "unknown")
{
    ## The message strings considered are the string constants subject to
    ## translation in gettext-family calls (see below for details).

    exprs <- parse(file = file, encoding = encoding, keep.source = TRUE)
    if(!length(exprs)) return(NULL)

    pd <- getParseData(exprs)

    ## Function for computing grandparent ids.
    parents <- pd$parent
    names(parents) <- pd$id
    gpids <- function(ids)
        parents[as.character(parents[as.character(ids)])]

    ind <- (pd$token == "SYMBOL_FUNCTION_CALL") &
        !is.na(match(pd$text,
                     c("warning", "stop",
                       "message", "packageStartupMessage",
                       "gettext", "gettextf", "ngettext")))

    funs <- pd$text[ind]

    ids <- gpids(pd$id[ind])
    calls <- getParseText(pd, ids)

    table <- pd[pd$token == "STR_CONST", ]
    pos <- match(gpids(table$id), ids)
    ind <- !is.na(pos)
    table <- split(table[ind, ], factor(pos[ind], seq_along(ids)))

    ## We have synopses
    ##   message(..., domain = NULL, appendLF = TRUE)
    ##   packageStartupMessage(..., domain = NULL, appendLF = TRUE)
    ##   warning(..., call. = TRUE, immediate. = FALSE, domain = NULL)
    ##   stop(..., call. = TRUE, domain = NULL)
    ##   gettext(..., domain = NULL)
    ##   ngettext(n, msg1, msg2, domain = NULL)
    ##   gettextf(fmt, ..., domain = NULL)
    ## For the first five, we simply take all unnamed strings.
    ## (Could make this more precise, of course.)
    ## For the latter two, we take the msg1/msg2 and fmt arguments,
    ## provided these are strings.

    ## <NOTE>
    ## Using domain = NA inhibits translation: perhaps it should
    ## optionally also inhibit spell checking?
    ## </NOTE>

    extract_message_strings <- function(fun, call, table) {
        ## Matching a call containing ... gives
        ##   Error in match.call(message, call) :
        ##   ... used in a situation where it doesn't exist
        ## so eliminate these.
        ## (Note that we also drop "..." strings.)
        call <- parse(text = call)[[1L]]
        call <- call[ as.character(call) != "..." ]
        mc <- as.list(match.call(get(fun, envir = .BaseNamespaceEnv),
                                 call))
        args <- if(fun == "gettextf")
            mc["fmt"]
        else if(fun == "ngettext")
            mc[c("msg1", "msg2")]
        else {
            if(!is.null(names(mc)))
                mc <- mc[!nzchar(names(mc))]
            mc[-1L]
        }
        strings <- as.character(args[vapply(args, is.character, TRUE)])
        ## Need to canonicalize to match string constants before and
        ## after parsing ...
        texts <- vapply(parse(text = table$text), as.character, "")
        pos <- which(!is.na(match(texts, strings)))
        cbind(table[pos, ], caller = rep.int(fun, length(pos)))
    }

    do.call(rbind,
            Map(extract_message_strings,
                as.list(funs), as.list(calls), table))
}

## For spell-checking the R R files.

aspell_R_R_files <-
function(which = NULL, dir = NULL,
         ignore = c("[ \t]'[^']*'[ \t[:punct:]]",
                    "[ \t][[:alnum:]_.]*\\(\\)[ \t[:punct:]]"),
         program = NULL, dictionaries = aspell_dictionaries_R)
{
    if(is.null(dir)) dir <- tools:::.R_top_srcdir_from_Rd()
    if(is.null(which))
        which <- tools:::.get_standard_package_names()$base

    files <-
        unlist(lapply(file.path(dir, "src", "library", which, "R"),
                      tools::list_files_with_type,
                      "code",
                      OS_subdirs = c("unix", "windows")),
               use.names = FALSE)

    program <- aspell_find_program(program)

    aspell(files,
           filter = list("R", ignore = ignore),
           control = aspell_control_R_Rd_files[[names(program)]],
           program = program,
           dictionaries = dictionaries)
}

## For spell-checking R files in a package.

aspell_package_R_files <-
function(dir, ignore = character(),
         control = list(), program = NULL, dictionaries = character())
{
    dir <- tools::file_path_as_absolute(dir)

    subdir <- file.path(dir, "R")
    files <- if(file_test("-d", subdir))
        tools::list_files_with_type(subdir,
                                    "code",
                                    OS_subdirs = c("unix", "windows"))
    else character()

    meta <- tools:::.get_package_metadata(dir, installed = FALSE)
    if(is.na(encoding <- meta["Encoding"]))
        encoding <- "unknown"

    defaults <- .aspell_package_defaults(dir, encoding)$R_files
    if(!is.null(defaults)) {
        if(!is.null(d <- defaults$ignore))
            ignore <- d
        if(!is.null(d <- defaults$control))
            control <- d
        if(!is.null(d <- defaults$program))
            program <- d
        if(!is.null(d <- defaults$dictionaries)) {
            dictionaries <-
                aspell_find_dictionaries(d, file.path(dir, ".aspell"))
        }
    }

    program <- aspell_find_program(program)

    aspell(files,
           filter = list("R", ignore = ignore),
           control = control,
           encoding = encoding,
           program = program,
           dictionaries = dictionaries)
}

## Spell-checking pot files.

## (Of course, directly analyzing the message strings would be more
## useful, but require writing appropriate text filters.)

## See also tools:::checkPoFile().

aspell_filter_db$pot <-
function (ifile, encoding = "unknown", ignore = character())
{
    lines <- readLines(ifile, encoding = encoding)

    ind <- grepl("^msgid[ \t]", lines)

    do_entry <- function(s) {
        out <- character(length(s))
        i <- 1L
        out[i] <- blank_out_regexp_matches(s[i], "^msgid[ \t]+\"")
        while(grepl("^\"", s[i <- i + 1L]))
            out[i] <- sub("^\"", " ", s[i])
        if(grepl("^msgid_plural[ \t]", s[i])) {
            out[i] <- blank_out_regexp_matches(s[i], "^msgid_plural[ \t]+\"")
            while(grepl("^\"", s[i <- i + 1L]))
                out[i] <- sub("^\"", " ", s[i])
        }
        out
    }

    entries <- split(lines, cumsum(ind))
    lines <- c(character(length(entries[[1L]])),
               as.character(do.call(c, lapply(entries[-1L], do_entry))))

    lines <- sub("\"[ \t]*$", " ", lines)

    ## <FIXME>
    ## Could replace backslash escapes for blanks and percent escapes by
    ## blanks, similar to what the R text filter does.
    ## </FIXME>

    for(re in ignore[nzchar(ignore)])
        lines <- blank_out_regexp_matches(lines, re)

    lines
}

## For spell-checking all pot files in a package.

aspell_package_pot_files <-
function(dir, ignore = character(),
         control = list(), program = NULL, dictionaries = character())
{
    dir <- tools::file_path_as_absolute(dir)
    subdir <- file.path(dir, "po")
    files <- if(file_test("-d", subdir))
        Sys.glob(file.path(subdir, "*.pot"))
    else character()

    meta <- tools:::.get_package_metadata(dir, installed = FALSE)
    if(is.na(encoding <- meta["Encoding"]))
        encoding <- "unknown"

    program <- aspell_find_program(program)

    aspell(files,
           filter = list("pot", ignore = ignore),
           control = control,
           encoding = encoding,
           program = program,
           dictionaries = dictionaries)
}

## For spell-checking the R C files.

aspell_R_C_files <-
function(which = NULL, dir = NULL,
         ignore = c("[ \t]'[[:alnum:]_.]*'[ \t[:punct:]]",
                    "[ \t][[:alnum:]_.]*\\(\\)[ \t[:punct:]]"),
         program = NULL, dictionaries = aspell_dictionaries_R)
{
    if(is.null(dir)) dir <- tools:::.R_top_srcdir_from_Rd()
    if(is.null(which))
        which <- tools:::.get_standard_package_names()$base
    if(!is.na(pos <- match("base", which)))
        which[pos] <- "R"

    files <- sprintf("%s.pot",
                     file.path(dir, "src", "library",
                               which, "po", which))
    files <- files[file_test("-f", files)]

    program <- aspell_find_program(program)

    aspell(files,
           filter = list("pot", ignore = ignore),
           control = aspell_control_R_Rd_files[[names(program)]],
           program = program,
           dictionaries = dictionaries)
}

## For spell-checking package C files.

aspell_package_C_files <-
function(dir, ignore = character(),
         control = list(), program = NULL, dictionaries = character())
{
    dir <- tools::file_path_as_absolute(dir)
    ## Assume that the package C message template file is shipped as
    ## 'po/PACKAGE.pot'.
    files <- file.path(dir, "po",
                       paste(basename(dir), "pot", collapse = "."))
    files <- files[file_test("-f", files)]

    meta <- tools:::.get_package_metadata(dir, installed = FALSE)
    if(is.na(encoding <- meta["Encoding"]))
        encoding <- "unknown"

    defaults <- .aspell_package_defaults(dir, encoding)$C_files
    if(!is.null(defaults)) {
        if(!is.null(d <- defaults$ignore))
            ignore <- d
        if(!is.null(d <- defaults$control))
            control <- d
        if(!is.null(d <- defaults$program))
            program <- d
        if(!is.null(d <- defaults$dictionaries)) {
            dictionaries <-
                aspell_find_dictionaries(d, file.path(dir, ".aspell"))
        }
    }

    program <- aspell_find_program(program)

    aspell(files,
           filter = list("pot", ignore = ignore),
           control = control,
           encoding = encoding,
           program = program,
           dictionaries = dictionaries)
}

## Spell-checking DCF files.

aspell_filter_db$dcf <-
function(ifile, encoding, keep = c("Title", "Description"),
         ignore = character())
{
    lines <- readLines(ifile, encoding = encoding)
    line_has_tags <- grepl("^[^[:blank:]][^:]*:", lines)
    tags <- sub(":.*", "", lines[line_has_tags])
    lines <- split(lines, cumsum(line_has_tags))
    ind <- is.na(match(tags, keep))
    lines[ind] <- lapply(lines[ind], function(s) rep.int("", length(s)))
    lines <- unlist(lines, use.names = FALSE)
    for(re in ignore[nzchar(ignore)])
        lines <- blank_out_regexp_matches(lines, re)
    lines
}

## For spell-checking package DESCRIPTION files.

aspell_package_description <-
function(dir, ignore = character(),
         control = list(), program = NULL, dictionaries = character())
{
    dir <- tools::file_path_as_absolute(dir)
    files <- file.path(dir, "DESCRIPTION")

    meta <- tools:::.get_package_metadata(dir, installed = FALSE)
    if(is.na(encoding <- meta["Encoding"]))
        encoding <- "unknown"

    program <- aspell_find_program(program)

    aspell(files,
           filter = list("dcf", ignore = ignore),
           control = control,
           encoding = encoding,
           program = program,
           dictionaries = dictionaries)
}

## For writing personal dictionaries:

aspell_write_personal_dictionary_file <-
function(x, out, language = "en", program = NULL)
{
    if(inherits(x, "aspell"))
        x <- sort(unique(x$Original))

    program <- aspell_find_program(program)
    if(is.na(program))
        stop("No suitable spell check program found.")

    ## <NOTE>
    ## Ispell and Hunspell take simple word lists as personal dictionary
    ## files, but Aspell requires a special format, see e.g.
    ## http://aspell.net/man-html/Format-of-the-Personal-and-Replacement-Dictionaries.html
    ## and one has to create these by hand, as
    ##   aspell --lang=en create personal ./foo "a b c"
    ## gives: Sorry "create/merge personal" is currently unimplemented.

    ## Encodings are a nightmare.
    ## Try to canonicalize to UTF-8 for Aspell (which allows recording
    ## the encoding in the personal dictionary).
    ## <FIXME>
    ## What should we do for Hunspell (which can handle UTF-8, but has
    ## no encoding information in the personal dictionary), or Ispell
    ## (which cannot handle UTF-8)?
    ## </FIXME>

    if(names(program) == "aspell") {
        header <- sprintf("personal_ws-1.1 %s %d UTF-8",
                          language, length(x))
        x <- enc2utf8(x)
    }
    else {
        header <- NULL
    }

    writeLines(c(header, x), out, useBytes = TRUE)
}

## For reading package defaults:

.aspell_package_defaults <-
function(dir, encoding = "unknown")
{
    dfile <- file.path(dir, ".aspell", "defaults.R")
    if(!file_test("-f", dfile))
        return(NULL)
    exprs <- parse(dfile, encoding = encoding)
    envir <- new.env()
    for(e in exprs) eval(e, envir)
    as.list(envir)
}

## Utilities.

blank_out_regexp_matches <-
function(s, re)
{
    m <- gregexpr(re, s)
    regmatches(s, m) <- Map(blanks, lapply(regmatches(s, m), nchar))
    s
}

blanks <-
function(n) {
    vapply(Map(rep.int, rep.int(" ", length(n)), n, USE.NAMES = FALSE),
           paste, "", collapse = "")
}

find_files_in_directories <-
function(basenames, dirnames)
{
    dirnames <- dirnames[file_test("-d", dirnames)]
    dirnames <- normalizePath(dirnames, "/")

    out <- character(length(basenames))
    pos <- seq_along(out)

    for(dir in dirnames) {
        paths <- file.path(dir, basenames[pos])
        ind <- file_test("-f", paths)
        out[pos[ind]] <- paths[ind]
        pos <- pos[!ind]
        if(!length(pos)) break
    }

    out
}
#  File src/library/utils/R/browseVignettes.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/



browseVignettes <- function(package = NULL, lib.loc = NULL, all = TRUE)
{
   
    vinfo <- tools:::getVignetteInfo(package, lib.loc, all)
    pkgs <- unique(vinfo[, "Package"])
    db <- lapply(pkgs, function(p) vinfo[vinfo[,"Package"] == p,,drop=FALSE])
    names(db) <- pkgs
    attr(db, "call") <- sys.call()
    attr(db, "footer") <-
        if (all) ""
        else sprintf(gettext("Use <code> %s </code> \n to list the vignettes in all <strong>available</strong> packages."),
                     "browseVignettes(all = TRUE)")
    class(db) <- "browseVignettes"
    return(db)
}

print.browseVignettes <- function(x, ...)
{
    if (length(x) == 0L) {
        message(gettextf("No vignettes found by %s",
                         paste(deparse(attr(x, "call")), collapse=" ")),
                domain = NA)
        return(invisible(x))
    }

    oneLink <- function(s) {
        if (length(s) == 0L) return(character(0L))
        title <- s[, "Title"]
        if (tools:::httpdPort > 0L)
            prefix <- sprintf("/library/%s/doc", pkg)
        else
            prefix <- sprintf("file://%s/doc", s[, "Dir"])
        src <- s[, "File"]
        pdf <- s[, "PDF"]
        rcode <- s[, "R"]
        pdfext <- sub("^.*\\.", "", pdf)
        sprintf("  <li>%s  -  \n    %s  \n    %s  \n    %s \n  </li>\n",
                title,
                ifelse(nzchar(pdf),
                       sprintf("<a href='%s/%s'>%s</a>&nbsp;",
                               prefix, pdf, toupper(pdfext)),
                       ""),
		sprintf("<a href='%s/%s'>source</a>&nbsp;", prefix, src),
		ifelse(nzchar(rcode),
                       sprintf("<a href='%s/%s'>R code</a>&nbsp;", prefix, rcode),
                       ""))
    }

    if (tools:::httpdPort == 0L)
        tools::startDynamicHelp()

    file <- tempfile("Rvig.", fileext=".html")
    sink(file)
    if (tools:::httpdPort > 0)
    	css_file <- "/doc/html/R.css"
    else
    	css_file <- file.path(R.home("doc"), "html", "R.css")
    cat(sprintf("<!DOCTYPE html PUBLIC '-//W3C//DTD HTML 4.01 Transitional//EN'>
<html>
<head>
<title>R Vignettes</title>
<meta http-equiv='Content-Type' content='text/html; charset=iso-8859-1'>
<link rel='stylesheet' type='text/css' href='%s'>
</head>
<body>\n", css_file))
    cat(sprintf("<h2>Vignettes found by <code><q>%s</q></code></h2>",
                paste(deparse(attr(x, "call")), collapse=" ")))
    cat("<div class=\"vignettes\">")
    for (pkg in names(x))
    {
        cat(sprintf("<h3>Vignettes in package <code>%s</code></h3>\n", pkg))
        cat("<ul>\n")
        links <- oneLink(x[[pkg]])
        cat(paste(links), collapse = "\n")
        cat("\n</ul>\n")
    }
    cat("</div>")
    cat(sprintf("<hr/><p>%s</p>", attr(x, "footer")))
    cat("</body></html>\n")
    sink()
    ## the first two don't work on Windows with browser=NULL.
    ## browseURL(URLencode(sprintf("file://%s", file)))
    ## browseURL(URLencode(file))
    if (tools:::httpdPort > 0L)
	browseURL(sprintf("http://127.0.0.1:%d/session/%s", tools:::httpdPort, basename(file)))
    else
    	browseURL(sprintf("file://%s", file))
    ## browseURL(file)
    invisible(x)
}
#  File src/library/utils/R/unix/bug.report.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

bug.report.info <- function()
    c("R Version:",
      paste0(" ", names(R.version), " = ", R.version),
      if (nzchar(Sys.getenv("R_GUI_APP_VERSION")))
          c("", "GUI:",
            paste0(" R-GUI ", Sys.getenv("R_GUI_APP_VERSION"),
                   " (", Sys.getenv("R_GUI_APP_REVISION"),")")),
      if (.Platform$OS.type == "windows") c("", win.version()),
      "",
      "Locale:", paste0(" ", Sys.getlocale()),
      "",
      "Search Path:",
      strwrap(paste(search(), collapse=", "), indent = 1, exdent = 1),
      "")

bug.report <- function(subject = "", address,
                       file = "R.bug.report", package = NULL, lib.loc = NULL,
                       ...)
{
    baseR <- function() {
        writeLines(c("  Bug reports on R and the base packages need to be submitted",
                     "  to the tracker at http://bugs.r-project.org/ .",
                     "",
                     "  We will now try to open that website in a browser"))
        flush.console()
        Sys.sleep(2)
        browseURL("https://bugs.r-project.org/bugzilla3/index.cgi")
    }

    findEmail <- function(x) {
        ## extract the part within the first < >: the rest may be invalid.
        x <- paste(x, collapse = " ") # could be multiple lines
        sub("[^<]*<([^>]+)>.*", "\\1", x)
    }
    if (is.null(package)) return(baseR())

    DESC <- packageDescription(package, lib.loc)
    if (!inherits(DESC, "packageDescription"))
        stop(gettextf("Package %s: DESCRIPTION file not found",
                      sQuote(package)), domain = NA)
    info <- paste0(c("Package", " Version", " Maintainer", " Built"),
		   ": ",
		   c(DESC$Package, DESC$Version, DESC$Maintainer, DESC$Built))
    info <- c(info, "", bug.report.info())
    if(identical(DESC$Priority, "base")) return(baseR())

    if (!is.null(DESC$BugReports)) {
        writeLines(info)
        cat("\nThis package has a bug submission web page, which we will now attempt\n",
            "to open.  The information above may be useful in your report. If the web\n",
            "page doesn't work, you should send email to the maintainer,\n",
            DESC$Maintainer, ".\n",
            sep = "")
        flush.console()
        Sys.sleep(2)
        browseURL(DESC$BugReports)
        return(invisible())
    }

    if (missing(address)) address <- findEmail(DESC$Maintainer)
    create.post(instructions = c("", "<<insert bug report here>>", rep("", 3)),
                description = "bug report",
                subject = subject, address = address,
                filename = file, info = info, ...)
}
#  File src/library/utils/R/capture.output.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

capture.output <- function(..., file=NULL, append=FALSE)
{
    args <- substitute(list(...))[-1L]

    rval <- NULL; closeit <- TRUE
    if (is.null(file))
        file <- textConnection("rval", "w", local = TRUE)
    else if (is.character(file))
        file <- file(file, if(append) "a" else "w")
    else if (inherits(file, "connection")) {
	if (!isOpen(file)) open(file, if(append) "a" else "w")
	else closeit <- FALSE
    } else
        stop("'file' must be NULL, a character string or a connection")

    sink(file)
    ## for error recovery: all output will be lost if file=NULL
    on.exit({sink(); if(closeit) close(file)})

    pf <- parent.frame()
    evalVis <- function(expr)
        withVisible(eval(expr, pf))

    for(i in seq_along(args)) {
        expr <- args[[i]]
        tmp <- switch(mode(expr),
                      "expression" = lapply(expr, evalVis),
                      "call" =, "name" =  list(evalVis(expr)),
                       stop("bad argument"))
        for(item in tmp)
            if (item$visible) print(item$value)
    }
    ## we need to close the text connection before returning 'rval'
    on.exit()
    sink()
    if(closeit) close(file)
    if(is.null(rval)) invisible(NULL) else rval
}
#  File src/library/utils/R/citation.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

## What a silly name ...
.is_not_nonempty_text <-
function(x)
    is.null(x) || any(is.na(x)) || all(grepl("^[[:space:]]*$", x))

person <-
function(given = NULL, family = NULL, middle = NULL,
         email = NULL, role = NULL, comment = NULL,
         first = NULL, last = NULL)
{
    ## Arrange all arguments in lists of equal length.
    args <- list(given = given, family = family, middle = middle,
                 email = email, role = role, comment = comment,
		 first = first, last = last)
    args <- lapply(args, .listify)
    args_length <- sapply(args, length)
    if(!all(args_length_ok <- args_length %in% c(1L, max(args_length))))
        warning(gettextf("Not all arguments are of the same length, the following need to be recycled: %s",
                         paste(names(args)[!args_length_ok],
                               collapse = ", ")),
                domain = NA)
    args <- lapply(args, function(x) rep(x, length.out = max(args_length)))

    ## <COMMENT Z>
    ## We could do this more elegantly, but let's just go through the
    ## list processing person by person.  I'm just recycling the old
    ## person() code for this.
    ## </COMMENT>
    person1 <-
    function(given = NULL, family = NULL, middle = NULL,
             email = NULL, role = NULL, comment = NULL,
             first = NULL, last = NULL)
    {
        if(!.is_not_nonempty_text(first)) {
            if(!.is_not_nonempty_text(given))
                stop(gettextf("Use either %s or %s/%s but not both.",
                              sQuote("given"),
                              sQuote("first"), sQuote("middle")),
                     domain = NA)
            ## <FIXME>
            ## Start warning eventually ... maybe use message() for now?
            message(gettextf("It is recommended to use %s instead of %s.",
                             sQuote("given"), sQuote("first")),
                    domain = NA)
            ## </FIXME>
            given <- first
        }
        if(!.is_not_nonempty_text(middle)) {
            ## <FIXME>
            ## Start warning eventually ... maybe use message() for now?
            message(gettextf("It is recommended to use %s instead of %s.",
                             sQuote("given"), sQuote("middle")),
                    domain = NA)
            ## </FIXME>
            given <- c(given, unlist(strsplit(middle, "[[:space:]]+")))
        }

        if(!.is_not_nonempty_text(last)) {
            if(!.is_not_nonempty_text(family))
                stop(gettextf("Use either %s or %s but not both.",
                              sQuote("family"), sQuote("last")),
                     domain = NA)
            ## <FIXME>
            ## Start warning eventually ... maybe use message() for now?
            message(gettextf("It is recommended to use %s instead of %s.",
                             sQuote("family"), sQuote("last")),
                    domain = NA)
            ## </FIXME>
            family <- last
        }

        ## Set all empty arguments to NULL.
        if(.is_not_nonempty_text(given)) given <- NULL
        if(.is_not_nonempty_text(family)) family <- NULL
        if(.is_not_nonempty_text(email)) email <- NULL
        if(.is_not_nonempty_text(role)) {
            if(!is.null(role))
                warning(sprintf(ngettext(length(role),
                                         "Invalid role specification: %s.",
                                         "Invalid role specifications: %s."),
                                paste(sQuote(role), collapse = ", ")),
                        domain = NA)
            role <- NULL
        }
        if(.is_not_nonempty_text(comment)) comment <- NULL

        ## <FIXME>
        ## Use something along the lines of
        ##   tools:::.valid_maintainer_field_regexp
        ## to validate given email addresses.
        ## </FIXME>

        if(length(role))
            role <- .canonicalize_person_role(role)

        rval <- list(given = given, family = family, role = role,
                     email = email, comment = comment)
        ## Canonicalize 0-length character arguments to NULL.
        if(any(ind <- (sapply(rval, length) == 0L)))
            rval[ind] <- vector("list", length = sum(ind))

        return(rval)
    }

    rval <-
        lapply(seq_along(args$given),
               function(i)
               with(args,
                    person1(given = given[[i]], family = family[[i]],
                            middle = middle[[i]], email = email[[i]],
                            role = role[[i]], comment = comment[[i]],
                            first = first[[i]], last = last[[i]])))
    class(rval) <- "person"

    ## <COMMENT Z>
    ## Should we check that for each person there is at least one
    ## non-NULL entry?
    ## </COMMENT>

    rval
}

.canonicalize_person_role <-
function(role)
{
    ## Be nice.  Given roles must either exactly match the role code,
    ## or be uniquely pmatchable modulo case against the role terms.
    pos <- which(is.na(match(role, MARC_relator_db$code)))
    if(length(pos)) {
        ind <- pmatch(tolower(role[pos]),
                      tolower(MARC_relator_db$name),
                      0L)
        role[pos[ind > 0L]] <- MARC_relator_db$code[ind]
        if(any(ind <- (ind == 0L))) {
            warning(sprintf(ngettext(length(pos[ind]),
                                     "Invalid role specification: %s.",
                                     "Invalid role specifications: %s."),
                            paste(sQuote(role[pos[ind]]), collapse = ", ")),
                    domain = NA)
            role <- role[-pos[ind]]
        }
    }
    role
}

`[[.person` <-
`[.person` <-
function(x, i)
{
    rval <- unclass(x)[i]
    class(rval) <- class(x)
    return(rval)
}

print.person <-
function(x, ...)
{
    x_char <- sapply(X = x, FUN = format, ...)
    print(x_char)
    invisible(x)
}

`$.person` <-
function(x, name)
{
    ## <COMMENT Z>
    ## extract internal list elements, return list if length > 1, vector
    ## otherwise (to mirror the behaviur of the input format for
    ## person())
    ## </COMMENT>
    name <- match.arg(name,
                      c("given", "family", "role", "email", "comment",
                        "first", "last", "middle")) # for now ...
    ## <COMMENT Z>
    ## Let's be nice and support first/middle/last for now.
    ## </COMMENT>
    if(name %in% c("first", "last", "middle")) {
        message(gettextf("It is recommended to use %s/%s instead of %s/%s/%s.",
                         sQuote("given"), sQuote("family"),
                         sQuote("first"), sQuote("middle"), sQuote("last")),
                domain = NA)
        oname <- name
	name <- switch(name,
	    "first" = "given",
	    "middle" = "given",
	    "last" = "family"
	)
    } else {
        oname <- name
    }

    rval <- lapply(unclass(x), function(p) p[[name]])

    if(oname == "first") rval <- lapply(rval, head, 1L)
    if(oname == "middle") {
        rval <- lapply(rval, tail, -1L)
        if(any(ind <- (sapply(rval, length) == 0L)))
            rval[ind] <- vector("list", length = sum(ind))
    }

    if(length(rval) == 1L) rval <- rval[[1L]]
    rval
}

`$<-.person` <-
function(x, name, value)
{
    name <- match.arg(name, c("given", "family", "role", "email", "comment"))
    x <- .listify(unclass(x))
    value <- rep(value, length.out = length(x))

    if(name == "role")
        value <- lapply(value, .canonicalize_person_role)

    for(i in seq_along(x)) {
        x[[i]][[name]] <- if(.is_not_nonempty_text(value[[i]]))
            NULL
        else as.character(value[[i]])
    }

    class(x) <- "person"
    x
}

c.person <-
function(..., recursive = FALSE)
{
    args <- list(...)
    if(!all(sapply(args, inherits, "person")))
        warning(gettextf("method is only applicable to %s objects",
                         sQuote("person")),
                domain = NA)
    args <- lapply(args, unclass)
    rval <- do.call("c", args)
    class(rval) <- "person"
    rval
}

as.person <-
function(x)
    UseMethod("as.person")

as.person.default <-
function(x)
{
    if(inherits(x, "person")) return(x)

    x <- as.character(x)

    ## Need to split the strings into individual person components.
    ## We used to split at ',' and 'and', but of course these could be
    ## contained in roles or comments as well.
    ## Hence, try the following.
    ## A. Replace all comment, role and email substrings by all-z
    ##    substrings of the same length.
    ## B. Tokenize the strings according to the split regexp matches in
    ##    the corresponding z-ified strings.
    ## C. Extract the persons from the thus obtained tokens.

    ## Create strings consisting of a given character c with given
    ## numbers n of characters.
    strings <- function(n, c = "z") {
        vapply(Map(rep.int, rep.int(c, length(n)), n,
                   USE.NAMES = FALSE),
               paste, "", collapse = "")
    }

    ## Replace matches of pattern in x by all-z substrings of the same
    ## length.
    zify <- function(pattern, x) {
        if(!length(x)) return(character())
        m <- gregexpr(pattern, x)
        regmatches(x, m) <-
            Map(strings, lapply(regmatches(x, m), nchar))
        x
    }

    ## Step A.
    y <- zify("\\([^)]*\\)", x)
    y <- zify("\\[[^]]*\\]", y)
    y <- zify("<[^>]*>", y)

    ## Step B.
    pattern <- "[[:space:]]?(,|,?[[:space:]]and)[[:space:]]+"
    x <- do.call("c",
                 regmatches(x, gregexpr(pattern, y), invert = TRUE))
    x <- x[!sapply(x, .is_not_nonempty_text)]

    ## Step C.
    as_person1 <- function(x) {
        comment <- if(grepl("\\(.*\\)", x))
            sub(".*\\(([^)]*)\\).*", "\\1", x)
        else NULL
        x <- sub("[[:space:]]*\\([^)]*\\)", "", x)
        email <- if(grepl("<.*>", x))
            sub(".*<([^>]*)>.*", "\\1", x)
        else NULL
        x <- sub("[[:space:]]*<[^>]*>", "", x)
        role <- if(grepl("\\[.*\\]", x))
            unlist(strsplit(gsub("[[:space:]]*", "",
                                 sub(".*\\[([^]]*)\\].*", "\\1", x)),
                            ",", fixed = TRUE))
        else NULL
        x <- sub("[[:space:]]*\\[[^)]*\\]", "", x)
        x <- unlist(strsplit(x, "[[:space:]]+"))
        z <- person(given = x[-length(x)], family = x[length(x)],
                    email = email, role = role, comment = comment)
        return(z)
    }

    as.list(do.call("c", lapply(x, as_person1)))
}

personList <-
function(...)
{
    z <- list(...)
    if(!all(sapply(z, inherits, "person")))
        stop(gettextf("all arguments must be of class %s",
                      dQuote("person")),
             domain = NA)
    do.call("c", z)
}

as.personList <-
function(x)
    UseMethod("as.personList")

as.personList.person <-
function(x)
    x

as.personList.default <-
function(x)
{
    if(inherits(x, "person")) return(x)
    do.call("c", lapply(x, as.person))
}

format.person <-
function(x,
         include = c("given", "family", "email", "role", "comment"),
         braces =
         list(given = "", family = "", email = c("<", ">"),
              role = c("[", "]"), comment = c("(", ")")),
         collapse =
         list(given = " ", family = " ", email = ", ",
              role = ", ", comment = ", "),
         ...
         )
{
    args <- c("given", "family", "email", "role", "comment")
    include <- sapply(include, match.arg, args)

    ## process defaults
    braces <- braces[args]
    collapse <- collapse[args]
    names(braces) <- names(collapse) <- args
    if(is.null(braces$given)) braces$given <- ""
    if(is.null(braces$family)) braces$family <- ""
    if(is.null(braces$email)) braces$email <- c("<", ">")
    if(is.null(braces$role)) braces$role <- c("[", "]")
    if(is.null(braces$comment)) braces$comment <- c("(", ")")
    braces <- lapply(braces, rep, length.out = 2L)
    if(is.null(collapse$given)) collapse$given <- " "
    if(is.null(collapse$family)) collapse$family <- " "
    if(is.null(collapse$email)) collapse$email <- ", "
    if(is.null(collapse$role)) collapse$role <- ", "
    if(is.null(collapse$comment)) collapse$comment <- ", "
    collapse <- lapply(collapse, rep, length.out = 1L)

    ## extract selected elements
    x <- lapply(unclass(x), "[", include)
    braces <- braces[include]
    collapse <- collapse[include]

    ## format 1 person
    format_person1 <- function(p) {
	rval <- lapply(seq_along(p), function(i) if(is.null(p[[i]])) NULL else
		       paste0(braces[[i]][1L], paste(p[[i]], collapse = collapse[[i]]),
			      braces[[i]][2L]))
	paste(do.call("c", rval), collapse = " ")
    }

    sapply(x, format_person1)
}

as.character.person <-
function(x, ...)
    format(x, ...)

toBibtex.person <-
function(object, ...)
    paste(format(object, include = c("given", "family")),
          collapse = " and ")

######################################################################

bibentry <-
function(bibtype, textVersion = NULL, header = NULL, footer = NULL, key = NULL,
         ...,
         other = list(), mheader = NULL, mfooter = NULL)
{
    BibTeX_names <- names(tools:::BibTeX_entry_field_db)

    args <- c(list(...), other)
    if(!length(args))
        return(structure(list(), class = "bibentry"))
    if(any(sapply(names(args), .is_not_nonempty_text)))
        stop("all fields have to be named")

    ## arrange all arguments in lists of equal length
    args <- c(list(bibtype = bibtype, textVersion = textVersion,
              header = header, footer = footer, key = key), list(...))
    args <- lapply(args, .listify)
    other <- lapply(other, .listify)
    max_length <- max(sapply(c(args, other), length))

    args_length <- sapply(args, length)
    if(!all(args_length_ok <- args_length %in% c(1L, max_length)))
        warning(gettextf("Not all arguments are of the same length, the following need to be recycled: %s",
                         paste(names(args)[!args_length_ok],
                               collapse = ", ")),
                domain = NA)
    args <- lapply(args, function(x) rep(x, length.out = max_length))

    other_length <- sapply(other, length)
    if(!all(other_length_ok <- other_length %in% c(1L, max_length)))
        warning(gettextf("Not all arguments are of the same length, the following need to be recycled: %s",
                         paste(names(other)[!other_length_ok],
                               collapse = ", ")),
                domain = NA)
    other <- lapply(other, function(x) rep(x, length.out = max_length))

    bibentry1 <-
    function(bibtype, textVersion, header = NULL, footer = NULL, key = NULL, ..., other = list())
    {
        ## process bibtype
	bibtype <- as.character(bibtype)
	stopifnot(length(bibtype) == 1L)
        pos <- match(tolower(bibtype), tolower(BibTeX_names))
	if(is.na(pos))
            stop(gettextf("%s has to be one of %s",
                          sQuote("bibtype"),
                          paste(BibTeX_names, collapse = ", ")),
                 domain = NA)
	bibtype <- BibTeX_names[pos]

        ## process fields
        rval <- c(list(...), other)
        rval <- rval[!sapply(rval, .is_not_nonempty_text)]
	fields <- tolower(names(rval))
        names(rval) <- fields
        attr(rval, "bibtype") <- bibtype

        ## check required fields
        .bibentry_check_bibentry1(rval)

        ## canonicalize
        pos <- fields %in% c("author", "editor")
	if(any(pos)) {
            for(i in which(pos)) rval[[i]] <- as.person(rval[[i]])
	}
	if(any(!pos)) {
            for(i in which(!pos)) rval[[i]] <- as.character(rval[[i]])
	}

        ## set attributes
        attr(rval, "key") <-
            if(is.null(key)) NULL else as.character(key)
        if(!is.null(textVersion))
            attr(rval, "textVersion") <- as.character(textVersion)
        if(!.is_not_nonempty_text(header))
            attr(rval, "header") <- paste(header, collapse = "\n")
        if(!.is_not_nonempty_text(footer))
            attr(rval, "footer") <- paste(footer, collapse = "\n")

        return(rval)
    }

    rval <- lapply(seq_along(args$bibtype),
                   function(i)
                   do.call("bibentry1",
                           c(lapply(args, "[[", i),
                             list(other = lapply(other, "[[", i)))))

    ## add main header/footer for overall bibentry vector
    if(!.is_not_nonempty_text(mheader))
        attr(rval, "mheader") <- paste(mheader, collapse = "\n")
    if(!.is_not_nonempty_text(mfooter))
        attr(rval, "mfooter") <- paste(mfooter, collapse = "\n")

    class(rval) <- "bibentry"
    rval
}

.bibentry_check_bibentry1 <-
function(x, force = FALSE)
{
    fields <- names(x)
    if(!force && !.is_not_nonempty_text(x$crossref)) return(NULL)
    bibtype <- attr(x, "bibtype")
    rfields <-
        strsplit(tools:::BibTeX_entry_field_db[[bibtype]], "|",
                 fixed = TRUE)
    if(length(rfields) > 0L) {
        ok <- sapply(rfields, function(f) any(f %in% fields))
        if(any(!ok))
            stop(sprintf(ngettext(sum(!ok),
                                  "A bibentry of bibtype %s has to specify the field: %s",
                                  "A bibentry of bibtype %s has to specify the fields: %s"),
                          sQuote(bibtype), paste(rfields[!ok], collapse = ", ")),
                 domain = NA)
    }
}

bibentry_attribute_names <-
    c("bibtype", "textVersion", "header", "footer", "key")
    
bibentry_list_attribute_names <-
    c("mheader", "mfooter")

.bibentry_get_key <-
function(x)
{
    if(!length(x)) return(character())
    keys <- lapply(unclass(x), attr, "key")
    keys[!vapply(keys, length, 0L)] <- ""
    unlist(keys)
}

`[[.bibentry` <-
`[.bibentry` <-
function(x, i, drop = TRUE)
{
    if(!length(x)) return(x)

    cl <- class(x)
    class(x) <- NULL
    ## For character subscripting, use keys if there are no names.
    ## Note that creating bibentries does not add the keys as names:
    ## assuming that both can independently be set, we would need to
    ## track whether names were auto-generated or not.
    ## (We could consider providing a names() getter which returns given
    ## names or keys as used for character subscripting, though).
    if(is.character(i) && is.null(names(x)))
        names(x) <- .bibentry_get_key(x)
    y <- x[i]
    if(!all(ok <- sapply(y, length) > 0L)) {
        warning("subscript out of bounds")
        y <- y[ok]
    }
    if(!drop)
        attributes(y) <- attributes(x)[bibentry_list_attribute_names]
    class(y) <- cl
    y
}

bibentry_format_styles <-
    c("text", "Bibtex", "citation", "html", "latex", "textVersion", "R")

.bibentry_match_format_style <-
function(style)
{
    ind <- pmatch(tolower(style), tolower(bibentry_format_styles),
                  nomatch = 0L)
    if(all(ind == 0L))
        stop(gettextf("%s should be one of %s",
                      sQuote("style"),
                      paste(dQuote(bibentry_format_styles),
                            collapse = ", ")),
             domain = NA)
    bibentry_format_styles[ind]
}

format.bibentry <-
function(x, style = "text", .bibstyle = NULL,
         citation.bibtex.max = getOption("citation.bibtex.max", 1),
         sort = FALSE, ...)
{
    style <- .bibentry_match_format_style(style)

    if(sort) x <- sort(x, .bibstyle = .bibstyle)
    x$.index <- as.list(seq_along(x))

    .format_bibentry_via_Rd <- function(f) {
        out <- file()
        saveopt <- tools::Rd2txt_options(width = getOption("width"))
        on.exit({tools::Rd2txt_options(saveopt); close(out)})
        sapply(.bibentry_expand_crossrefs(x),
               function(y) {
                   rd <- tools::toRd(y, style = .bibstyle)
                   con <- textConnection(rd)
                   on.exit(close(con))
                   f(con, fragment = TRUE, out = out, ...)
                   paste(readLines(out), collapse = "\n")
               })
    }

    .format_bibentry_as_citation <- function(x) {
        bibtex <- length(x) <= citation.bibtex.max

        c(paste(strwrap(attr(x, "mheader")), collapse = "\n"),
          unlist(lapply(x, function(y) {
              paste(c(if(!is.null(y$header))
                      c(strwrap(y$header), ""),
                      if(!is.null(y$textVersion)) {
                          strwrap(y$textVersion, prefix = "  ")
                      } else {
                          format(y)
                      },
                      if(bibtex) {
                          c(gettext("\nA BibTeX entry for LaTeX users is\n"),
			    paste0("  ", unclass(toBibtex(y))))
                      },
                      if(!is.null(y$footer))
                      c("", strwrap(y$footer))),
                    collapse = "\n")
          })),
          paste(strwrap(attr(x, "mfooter")), collapse = "\n")
          )
    }

    out <-
        switch(style,
               "text" = .format_bibentry_via_Rd(tools::Rd2txt),
               "html" = .format_bibentry_via_Rd(tools::Rd2HTML),
               "latex" = .format_bibentry_via_Rd(tools::Rd2latex),
               "Bibtex" = {
                   unlist(lapply(x,
                                 function(y)
                                 paste(toBibtex(y), collapse = "\n")))
               },
               "textVersion" = {
                   out <- lapply(unclass(x), attr, "textVersion")
                   out[!sapply(out, length)] <- ""
                   unlist(out)
               },
               "citation" = .format_bibentry_as_citation(x),
               "R" = .format_bibentry_as_R_code(x, ...)
               )
    as.character(out)
}

.bibentry_expand_crossrefs <-
function(x, more = list())
{
    y <- if(length(more))
        do.call(c, c(list(x), more))
    else
        x

    x <- unclass(x)
    y <- unclass(y)

    crossrefs <- lapply(x, `[[`, "crossref")
    pc <- which(vapply(crossrefs, length, 0L) > 0L)

    if(length(pc)) {
        pk <- match(unlist(crossrefs[pc]), .bibentry_get_key(y))
        ## If an entry has a crossref we cannot resolve it might still
        ## be complete: we could warn about the bad crossref ...
        ok <- !is.na(pk)
        ## Merge entries: note that InCollection and InProceedings need
        ## to remap title to booktitle as needed.
        x[pc[ok]] <-
            Map(function(u, v) {
                add <- setdiff(names(v), names(u))
                u[add] <- v[add]
                if(!is.na(match(tolower(attr(u, "bibtype")),
                                c("incollection", "inproceedings"))) &&
                   is.null(u$booktitle))
                    u$booktitle <- v$title
                u
            },
                x[pc[ok]],
                y[pk[ok]])
        ## Now check entries with crossrefs for completeness.
        ## Ignore bad entries with a warning.
        status <- lapply(x[pc],
                         function(e)
                         tryCatch(.bibentry_check_bibentry1(e, TRUE),
                                  error = identity))
        bad <- which(sapply(status, inherits, "error"))
        if(length(bad)) {
            for(b in bad) {
                warning(gettextf("Dropping invalid entry %d:\n%s",
                                 pc[b],
                                 conditionMessage(status[[b]])))
            }
            x[pc[bad]] <- NULL
        }
    }

    class(x) <- "bibentry"
    x
}

print.bibentry <-
function(x, style = "text", .bibstyle = NULL, ...)
{
    style <- .bibentry_match_format_style(style)

    if(style == "R") {
	writeLines(format(x, "R", collapse = TRUE, ...))
    } else if(length(x)) {
	y <- format(x, style, .bibstyle, ...)
        if(style == "citation") {
            ## Printing in citation style does extra headers/footers
            ## (which however may be empty), so it is handled
            ## differently.
            n <- length(y)
            if(nzchar(header <- y[1L]))
                header <- c("", header, "")
            if(nzchar(footer <- y[n]))
                footer <- c("", footer, "")
            writeLines(c(header,
                         paste(y[-c(1L, n)], collapse = "\n\n"),
                         footer))
        } else {
            writeLines(paste(y, collapse = "\n\n"))
        }
    }

    invisible(x)
}

## Not vectorized for now: see ?regmatches for a vectorized version.
.blanks <-
function(n)
    paste(rep.int(" ", n), collapse = "")

.format_call_RR <-
function(cname, cargs)
{
    ## Format call with ragged right argument list (one arg per line).
    cargs <- as.list(cargs)
    n <- length(cargs)
    lens <- sapply(cargs, length)
    sums <- cumsum(lens)
    starters <- c(sprintf("%s(", cname),
                  rep.int(.blanks(nchar(cname) + 1L), sums[n] - 1L))
    trailers <- c(rep.int("", sums[n] - 1L), ")")
    trailers[sums[-n]] <- ","
    sprintf("%s%s%s", starters, unlist(cargs), trailers)
}

.format_bibentry_as_R_code <-
function(x, collapse = FALSE)
{
    if(!length(x)) return("bibentry()")

    x$.index <- NULL

    ## There are two subleties for constructing R calls giving a given
    ## bibentry object.
    ## * There can be mheader and mfooter entries.
    ##   If there are, we put them into the first bibentry.
    ## * There could be field names which clash with the names of the
    ##   bibentry() formals: these would need to be put as a list into
    ##   the 'other' formal.

    ## The following make it into the attributes of an entry.
    anames <- bibentry_attribute_names
    ## The following make it into the attributes of the object.
    manames <- c("mheader", "mfooter")

    ## Format a single element (person or string, at least for now).
    f <- function(e) {
        if(inherits(e, "person"))
            .format_person_as_R_code(e)
        else
            deparse(e)
    }

    g <- function(u, v) {
        prefix <- sprintf("%s = ", u)
        n <- length(v)
        if(n > 1L)
            prefix <- c(prefix,
                        rep.int(.blanks(nchar(prefix)), n - 1L))
        sprintf("%s%s", prefix, v)
    }

    s <- lapply(unclass(x),
                function(e) {
                    a <- Filter(length, attributes(e)[anames])
                    e <- e[!sapply(e, is.null)]
                    ind <- !is.na(match(names(e),
                                       c(anames, manames, "other")))
                    if(any(ind)) {
                        other <- paste(names(e[ind]),
                                       sapply(e[ind], f),
                                       sep = " = ")

                        other <- Map(g,
                                     names(e[ind]),
                                     sapply(e[ind], f))
                        other <- .format_call_RR("list", other)
                        e <- e[!ind]
                    } else {
                        other <- NULL
                    }
                    c(Map(g, names(a), sapply(a, deparse)),
                      Map(g, names(e), sapply(e, f)),
                      if(length(other)) list(g("other", other)))

                })

    if(!is.null(mheader <- attr(x, "mheader")))
        s[[1L]] <- c(s[[1L]],
                     paste("mheader = ", deparse(mheader)))
    if(!is.null(mfooter <- attr(x, "mfooter")))
        s[[1L]] <- c(s[[1L]],
                     paste("mfooter = ", deparse(mfooter)))

    s <- Map(.format_call_RR, "bibentry", s)
    if(collapse && (length(s) > 1L))
        paste(.format_call_RR("c", s), collapse = "\n")
    else
        unlist(lapply(s, paste, collapse = "\n"), use.names = FALSE)

}

.format_person_as_R_code <-
function(x)
{
    s <- lapply(unclass(x),
                function(e) {
                    e <- e[!sapply(e, is.null)]
                    cargs <-
                        sprintf("%s = %s", names(e), sapply(e, deparse))
                    .format_call_RR("person", cargs)
                })
    if(length(s) > 1L)
        .format_call_RR("c", s)
    else
        unlist(s, use.names = FALSE)
}

`$.bibentry` <-
function(x, name)
{
    if(!length(x)) return(NULL)

    ## <COMMENT Z>
    ## Extract internal list elements, return list if length > 1, vector
    ## otherwise (to mirror the behaviour of the input format for
    ## bibentry())
    ## </COMMENT>
    is_attribute <- name %in% bibentry_attribute_names
    rval <- if(is_attribute) lapply(unclass(x), attr, name)
        else lapply(unclass(x), "[[", name)
    if(length(rval) == 1L) rval <- rval[[1L]]
    rval
}

`$<-.bibentry` <-
function(x, name, value)
{
    is_attribute <- name %in% bibentry_attribute_names

    x <- unclass(x)
    name <- tolower(name)

    ## recycle value
    value <- rep(.listify(value), length.out = length(x))

    ## check bibtype
    if(name == "bibtype") {
        stopifnot(all(sapply(value, length) == 1L))
        BibTeX_names <- names(tools:::BibTeX_entry_field_db)
        value <- unlist(value)
        pos <- match(tolower(value), tolower(BibTeX_names))
        if(any(is.na(pos)))
            stop(gettextf("%s has to be one of %s",
                          sQuote("bibtype"),
                          paste(BibTeX_names, collapse = ", ")),
                 domain = NA)
        value <- as.list(BibTeX_names[pos])
    }

    ## replace all values
    for(i in seq_along(x)) {
        if(is_attribute) {
	    attr(x[[i]], name) <-
                if(is.null(value[[i]])) NULL else paste(value[[i]])
	} else {
	    x[[i]][[name]] <-
                if(is.null(value[[i]])) NULL else {
                    if(name %in% c("author", "editor"))
                        as.person(value[[i]])
                    else paste(value[[i]])
                }
        }
    }

    ## check whether all elements still have their required fields
    for(i in seq_along(x)) .bibentry_check_bibentry1(x[[i]])

    class(x) <- "bibentry"
    x
}

c.bibentry <-
function(..., recursive = FALSE)
{
    args <- list(...)
    if(!all(sapply(args, inherits, "bibentry")))
        warning(gettextf("method is only applicable to %s objects",
                         sQuote("bibentry")),
                domain = NA)
    args <- lapply(args, unclass)
    rval <- do.call("c", args)
    class(rval) <- "bibentry"
    rval
}

toBibtex.bibentry <-
function(object, ...)
{
    format_author <- function(author) paste(sapply(author, function(p) {
	fnms <- p$family
	only_given_or_family <- is.null(fnms) || is.null(p$given)
	fbrc <- if(length(fnms) > 1L ||
                   any(grepl("[[:space:]]", fnms)) ||
                   only_given_or_family) c("{", "}") else ""
	gbrc <- if(only_given_or_family) c("{", "}") else ""
        format(p, include = c("given", "family"),
               braces = list(given = gbrc, family = fbrc))
    }), collapse = " and ")

    format_bibentry1 <- function(object) {
	object <- unclass(object)[[1L]]
        rval <- paste0("@", attr(object, "bibtype"), "{", attr(object, "key"), ",")
        if("author" %in% names(object))
            object$author <- format_author(object$author)
        if("editor" %in% names(object))
            object$editor <- format_author(object$editor)

        rval <- c(rval,
                  sapply(names(object), function (n)
                         paste0("  ", n, " = {", object[[n]], "},")),
                  "}", "")
        return(rval)
    }

    if(length(object)) {
        object$.index <- NULL
        rval <- head(unlist(lapply(object, format_bibentry1)), -1L)
    } else
        rval <- character()
    class(rval) <- "Bibtex"
    rval
}

sort.bibentry <-
function(x, decreasing = FALSE, .bibstyle = NULL, drop = FALSE, ...)
{
    x[order(tools::bibstyle(.bibstyle)$sortKeys(x),
            decreasing = decreasing),
      drop = drop]
}

rep.bibentry <-
function(x, ...)
{
    y <- NextMethod("rep")
    class(y) <- class(x)
    y
}

unique.bibentry <-
function(x, ...)
{
    y <- NextMethod("unique")
    class(y) <- class(x)
    y
}

######################################################################

citEntry <-
function(entry, textVersion, header = NULL, footer = NULL, ...)
    bibentry(bibtype = entry, textVersion = textVersion,
             header = header, footer = footer, ...)

citHeader <-
function(...)
{
    rval <- paste(...)
    class(rval) <- "citationHeader"
    rval
}

citFooter <-
function(...)
{
    rval <- paste(...)
    class(rval) <- "citationFooter"
    rval
}

readCitationFile <-
function(file, meta = NULL)
{
    exprs <- tools:::.parse_CITATION_file(file, meta$Encoding)

    rval <- list()
    mheader <- NULL
    mfooter <- NULL
    k <- 0L
    envir <- new.env(hash = TRUE)
    ## Make the package metadata available to the citation entries.
    assign("meta", meta, envir = envir)

    for(expr in exprs) {
        x <- eval(expr, envir = envir)
        if(inherits(x, "bibentry"))
            rval <- c(rval, list(x))
        else if(identical(class(x), "citationHeader"))
            mheader <- c(mheader, x)
        else if(identical(class(x), "citationFooter"))
            mfooter <- c(mfooter, x)
    }

    rval <- if(length(rval) == 1L)
        rval[[1L]]
    else
        do.call("c", rval)
    if(!.is_not_nonempty_text(mheader))
        attr(rval, "mheader") <- paste(mheader, collapse = "\n")
    if(!.is_not_nonempty_text(mfooter))
        attr(rval, "mfooter") <- paste(mfooter, collapse = "\n")

    .citation(rval)
}

######################################################################

citation <-
function(package = "base", lib.loc = NULL, auto = NULL)
{
    ## Allow citation(auto = meta) in CITATION files to include
    ## auto-generated package citation.
    if(inherits(auto, "packageDescription")) {
        auto_was_meta <- TRUE
        meta <- auto
        package <- meta$Package
    } else {
        auto_was_meta <- FALSE
        dir <- system.file(package = package, lib.loc = lib.loc)
        if(dir == "")
            stop(gettextf("package %s not found", sQuote(package)),
                 domain = NA)
        meta <- packageDescription(pkg = package,
                                   lib.loc = dirname(dir))
        ## if(is.null(auto)): Use default auto-citation if no CITATION
        ## available.
        citfile <- file.path(dir, "CITATION")
        if(is.null(auto)) auto <- !file_test("-f", citfile)
        ## if CITATION is available
        if(!auto) {
            return(readCitationFile(citfile, meta))
        } else if(package == "base") {
            ## Avoid infinite recursion for broken installation.
            stop("broken installation, no CITATION file in the base package.")
        }
    }

    ## Auto-generate citation info.

    ## Base packages without a CITATION file use the base citation.
    if((!is.null(meta$Priority)) && (meta$Priority == "base")) {
    	cit <- citation("base", auto = FALSE)
    	attr(cit, "mheader")[1L] <-
	    paste0("The ", sQuote(package), " package is part of R.  ",
		   attr(cit, "mheader")[1L])
        return(.citation(cit))
    }

    year <- sub("-.*", "", meta$`Date/Publication`)
    if(!length(year)) {
        year <- sub(".*((19|20)[[:digit:]]{2}).*", "\\1", meta$Date,
                    perl = TRUE) # may not be needed, but safer
        if(is.null(meta$Date)){
            warning(gettextf("no date field in DESCRIPTION file of package %s",
                             sQuote(package)),
                    domain = NA)
        }
        else if(!length(year)) {
            warning(gettextf("could not determine year for %s from package DESCRIPTION file",
                             sQuote(package)),
                    domain = NA)
        }
    }

    author <- meta$`Authors@R`
    ## <FIXME>
    ## Older versions took persons with no roles as "implied" authors.
    ## So for now check whether Authors@R gives any authors; if not fall
    ## back to the plain text Author field.
    if(length(author)) {
        author <- .read_authors_at_R_field(author)
        ## We only want those with author roles.
        author <- Filter(.person_has_author_role, author)
    }
    if(length(author)) {
        has_authors_at_R_field <- TRUE
    } else {
        has_authors_at_R_field <- FALSE
        author <- as.personList(meta$Author)
    }
    ## </FIXME>

    z <- list(title = paste0(package, ": ", meta$Title),
              author = author,
              year = year,
              note = paste("R package version", meta$Version)
              )

    z$url <- if(identical(meta$Repository, "CRAN"))
        sprintf("http://CRAN.R-project.org/package=%s", package)
    else
        meta$URL

    if(identical(meta$Repository, "R-Forge")) {
        z$url <- if(!is.null(rfp <- meta$"Repository/R-Forge/Project"))
            sprintf("http://R-Forge.R-project.org/projects/%s/", rfp)
        else
            "http://R-Forge.R-project.org/"
        if(!is.null(rfr <- meta$"Repository/R-Forge/Revision"))
            z$note <- paste(z$note, rfr, sep = "/r")
    }

    header <- if(!auto_was_meta) {
        gettextf("To cite package %s in publications use:",
                 sQuote(package))
    } else NULL


    ## No auto-generation message for auto was meta so that maintainers
    ## can safely use citation(auto = meta) in their CITATION without
    ## getting notified about possible needs for editing.
    footer <- if(!has_authors_at_R_field && !auto_was_meta) {
        gettextf("ATTENTION: This citation information has been auto-generated from the package DESCRIPTION file and may need manual editing, see %s.",
                 sQuote("help(\"citation\")"))
    } else NULL

    author <- format(z$author, include = c("given", "family"))
    if(length(author) > 1L)
        author <- paste(paste(head(author, -1L), collapse = ", "),
                        tail(author, 1L), sep = " and ")

    rval <- bibentry(bibtype = "Manual",
                     textVersion =
                     paste0(author, " (", z$year, "). ", z$title, ". ",
                            z$note, ". ", z$url),
                     header = header,
                     footer = footer,
                     other = z
                     )
    .citation(rval)
}

.citation <-
function(x)
{
    class(x) <- c("citation", "bibentry")
    x
}

.read_authors_at_R_field <-
function(x)
{
    out <- eval(parse(text = x))

    ## Let's by nice ...
    ## Alternatively, we could throw an error.
    if(!inherits(out, "person"))
        out <- do.call("c", lapply(x, as.person))

    out
}

.person_has_author_role <-
function(x)
{
    ## <NOTE>
    ## Earlier versions used
    ##    is.null(r <- x$role) || "aut" %in% r
    ## using author roles by default.
    ## </NOTE>
    "aut" %in% x$role
}

print.citation <-
function(x, style = "citation", ...)
{
    NextMethod("print", x, style = style, ...)
    invisible(x)
}

as.bibentry <-
function(x)
    UseMethod("as.bibentry")

as.bibentry.bibentry <- identity

as.bibentry.citation <-
function(x)
{
    class(x) <- "bibentry"
    x
}

.listify <-
function(x)
    if(inherits(x, "list")) x else list(x)

.format_person_for_plain_author_spec <-
function(x) {
    ## Names first.
    out <- format(x, include = c("given", "family"))
    ## Only show roles recommended for usage with R.
    role <- x$role
    if(!length(role)) return("")
    role <- role[role %in% MARC_relator_db_codes_used_with_R]
    if(!length(role)) return("")
    out <- sprintf("%s [%s]", out, paste(role, collapse = ", "))
    if(!is.null(comment <- x$comment))
        out <- sprintf("%s (%s)", out,
                       paste(comment, collapse = "\n"))
    out
}

## NB: because of the use of strwrap(), this always outputs
## in the current locale even if the input has a marked encoding.
.format_authors_at_R_field_for_author <-
function(x)
{
    if(is.character(x))
        x <- .read_authors_at_R_field(x)
    header <- attr(x, "header")
    footer <- attr(x, "footer")
    x <- sapply(x, .format_person_for_plain_author_spec)
    ## Drop persons with irrelevant roles.
    x <- x[x != ""]
    ## And format.
    if(!length(x)) return("")
    ## We need to ensure that the first line has no indentation, whereas
    ## all subsequent lines are indented (as .write_description avoids
    ## folding for Author fields).  We use a common indentation of 2,
    ## with an extra indentation of 2 within single author descriptions.
    out <- paste(lapply(strwrap(x, indent = 0L, exdent = 4L,
                                simplify = FALSE),
                        paste, collapse = "\n"),
                 collapse = ",\n  ")
    if(!is.null(header)) {
        header <- paste(strwrap(header, indent = 0L, exdent = 2L),
                        collapse = "\n")
        out <- paste(header, out, sep = "\n  ")
    }
    if(!is.null(footer)) {
        footer <- paste(strwrap(footer, indent = 2L, exdent = 2L),
                        collapse = "\n")
        out <- paste(out, footer, sep = ".\n")
    }
    out
}

## preserves encoding if marked.
.format_authors_at_R_field_for_maintainer <-
function(x)
{
    if(is.character(x))
        x <- .read_authors_at_R_field(x)
    ## Maintainers need cre roles and email addresses.
    x <- Filter(function(e)
                !is.null(e$email) && ("cre" %in% e$role),
                x)
    ## If this leaves nothing ...
    if(!length(x)) return("")
    paste(format(x, include = c("given", "family", "email")),
          collapse = ",\n  ")
}

# Cite using the default style (which is usually citeNatbib)

cite <-
function(keys, bib, ...)
{
    fn <- tools::bibstyle()$cite
    if (is.null(fn))
    	fn <- citeNatbib
    fn(keys, bib, ...)
}

# Cite using natbib-like options.  A bibstyle would normally
# choose some of these options and just have a cite(keys, bib, previous)
# function within it.

citeNatbib <-
local({
    cited <- c()

    function(keys, bib, textual = FALSE, before = NULL, after = NULL,
             mode = c("authoryear", "numbers", "super"),
             abbreviate = TRUE, longnamesfirst = TRUE,
             bibpunct = c("(", ")", ";", "a", "", ","),
             previous) {

	shortName <- function(person) {
	    if (length(person$family))
		paste(tools:::cleanupLatex(person$family), collapse = " ")
	    else
		paste(tools:::cleanupLatex(person$given), collapse = " ")
	}

	authorList <- function(paper)
	    names <- sapply(paper$author, shortName)

	if (!missing(previous))
	    cited <<- previous

	if (!missing(mode))
	    mode <- match.arg(mode)
	else
	    mode <- switch(bibpunct[4L],
	    	n = "numbers",
	    	s = "super",
	    	"authoryear")
        numeric <- mode %in% c('numbers', 'super')

	if (numeric)
	    bib <- sort(bib)

	keys <- unlist(strsplit(keys, " *, *"))
	if (!length(keys)) return("")

        n <- length(keys)
	first <- !(keys %in% cited)
	cited <<- unique(c(cited, keys))

	bibkeys <- unlist(bib$key)
	# Use year to hold numeric entry; makes things
	# simpler below
	year <- match(keys, bibkeys)
	papers <- bib[year]

        if (textual || !numeric) {
	    auth <- character(n)
	    if (!numeric)
	    	year <- unlist(papers$year)
	    authorLists <- lapply(papers, authorList)
	    lastAuthors <- NULL
	    for (i in seq_along(keys)) {
		authors <- authorLists[[i]]
		if (identical(lastAuthors, authors))
		    auth[i] <- ""
		else {
		    if (length(authors) > 1L)
			authors[length(authors)] <- paste("and", authors[length(authors)])
		    if (length(authors) > 2L) {
			if (!abbreviate || (first[i] && longnamesfirst))
			    auth[i] <- paste(authors, collapse=", ")
			else
			    auth[i] <- paste(authors[1L], "et al.")
		    } else
			auth[i] <- paste(authors, collapse=" ")
            	}
            	lastAuthors <- authors
            }
            suppressauth <- which(!nzchar(auth))
            if (length(suppressauth)) {
                for (i in suppressauth)
                    year[i - 1L] <-
                        paste0(year[i - 1L], bibpunct[6L], " ", year[i])
                auth <- auth[-suppressauth]
                year <- year[-suppressauth]
            }
        }
        if (!is.null(before))
            before <- paste0(before, " ")
        if (!is.null(after))
            after <- paste0(" ", after)
        if (textual) {
            result <- paste0(bibpunct[1L], before, year, after, bibpunct[2L])
            if (mode == "super")
            	result <- paste0(auth, "^{", result, "}")
            else
            	result <- paste0(auth, " ", result)
            result <- paste(result, collapse = paste0(bibpunct[3L], " "))
        } else if (numeric) {
            result <- paste(year, collapse=paste0(bibpunct[3L], " "))
            result <- paste0(bibpunct[1L], before, result, after, bibpunct[2L])
            if (mode == "super")
            	result <- paste0("^{", result, "}")
        } else {
            result <- paste0(auth, bibpunct[5L], " ", year)
            result <- paste(result, collapse = paste0(bibpunct[3L], " "))
            result <- paste0(bibpunct[1L], before, result, after, bibpunct[2L])
        }
        result
    }
})
#  File src/library/utils/R/combn.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

combn <- function(x, m, FUN = NULL, simplify = TRUE, ...)
{
    ## DATE WRITTEN: 14 April 1994	    LAST REVISED:  10 July 1995
    ## AUTHOR:	Scott Chasalow
    ##
    ## DESCRIPTION:
    ##	Generate all combinations of the elements of x taken m at a time.
    ##	If x is a positive integer,  returns all combinations
    ##	of the elements of seq(x) taken m at a time.
    ##	If argument "FUN" is not null,	applies a function given
    ##	by the argument to each point.	If simplify is FALSE,  returns
    ##	a list; else returns a vector or an array.  "..." are passed
    ##	unchanged to function given by argument FUN,  if any.

    ##S : Change if (simplify = TRUE) return an array/matrix {not a 'vector'}
    stopifnot(length(m) == 1L, is.numeric(m))
    if(m < 0) stop("m < 0", domain = NA)
    if(is.numeric(x) && length(x) == 1L && x > 0 && trunc(x) == x)
	x <- seq_len(x)
    n <- length(x)
    if(n < m) stop("n < m", domain = NA)
    m <- as.integer(m)
    e <- 0
    h <- m
    a <- seq_len(m)
    nofun <- is.null(FUN)
    if(!nofun && !is.function(FUN))
	stop("'FUN' must be a function or NULL")
    # first result : what kind, what length,.. ?
    len.r <- length(r <- if(nofun) x[a] else FUN(x[a], ...))
    count <- as.integer(round(choose(n, m))) # >= 1
    if(simplify) {
	dim.use <-
	    if(nofun)
		c(m, count) # matrix also when count = 1
	    else {
		d <- dim(r)
		if(length(d) > 1L)
		    c(d, count)
		else if(len.r > 1L)
		    c(len.r, count)
		else # MM: *still* a matrix - a la "drop = FALSE"
		    c(d, count)
	    } ## NULL in all 'else' cases
##S	use.arr <- !is.null(dim.use)
    }
##S	else use.arr <- FALSE

    if(simplify) { # use atomic vector/array instead of list
##S	if(use.arr)
	    out <- matrix(r, nrow = len.r, ncol = count) # matrix for now
##S	else {
##S	    if(count > 1) {
##S		out <- vector(storage.mode(r), len.r * count)
##S		out[1L] <- r
##S	    }
##S	    else out <- r
##S	}
    }
    else {
	out <- vector("list", count)
	out[[1L]] <- r
    }

    if(m > 0) {
	i <- 2L
	nmmp1 <- n - m + 1L	 # using 1L to keep integer arithmetic
	while(a[1L] != nmmp1) {
	    if(e < n - h) {
		h <- 1L
		e <- a[m]
		j <- 1L
	    }
	    else {
		e <- a[m - h]
		h <- h + 1L
		j <- 1L:h
	    }
	    a[m - h + j] <- e + j
	    r <- if(nofun) x[a] else FUN(x[a], ...)
	    if(simplify) ##S if(use.arr)
		out[, i] <- r else out[[i]] <- r
	    i <- i + 1L
	}
    }
    if(simplify) ##S if(use.arr)
	array(out, dim.use) else out
}
#  File src/library/utils/R/completion.R
#  Part of the R package, http://www.R-project.org
#
# Copyright (C) 2006  Deepayan Sarkar
#               2006-2012  The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/



### Note: By default, we try not to do things that might be slow due
### to network latency (think NFS).  For example, retrieving a list of
### available packages is potentially slow, and is thus disabled
### initially.


### Status: I'm mostly happy with things.  The only obvious
### improvement I can think of is figuring out when we are in
### continuation mode (R prompt == "+") and make use of previous lines
### in that case.  I haven't found a way to do that.



### Note: sprintf seems faster than paste based on naive benchmarking:

## > system.time(for (i in 1L:100000L) sprintf("foo%sbar%d", letters, 1L:26L) )
##            user          system           total   user.children system.children
##           4.796           0.088           4.887           0.000           0.000
## > system.time(for (i in 1L:100000L) paste("foo", letters, "bar", 1L:26L) )
##            user          system           total   user.children system.children
##           8.300           0.028           8.336           0.000           0.000

### so will change all pastes to sprintf.  However, we need to be
### careful because 0 length components in sprintf will cause errors.


## generic and built-in methods to generate completion after $

.DollarNames <- function(x, pattern)
    UseMethod(".DollarNames")

.DollarNames.default <- function(x, pattern = "") {
    if (is.atomic(x) || is.symbol(x)) character()
    else grep(pattern, names(x), value = TRUE)
}

.DollarNames.list <- function(x, pattern = "") {
    grep(pattern, names(x), value = TRUE)
}

.DollarNames.environment <- function(x, pattern = "") {
    ls(x, all.names = TRUE, pattern = pattern)
}

## if (is.environment(object))
## {
##     ls(object,
##        all.names = TRUE,
##        pattern = sprintf("^%s", makeRegexpSafe(suffix)))
## }
## else
## {
##     grep(sprintf("^%s", makeRegexpSafe(suffix)),
##          names(object), value = TRUE)
## }




## modifies settings:

rc.settings <- function(ops, ns, args, func, ipck, S3, data, help, argdb, quotes, files)
{
    if (length(match.call()) == 1) return(unlist(.CompletionEnv[["settings"]]))
    checkAndChange <- function(what, value)
    {
        if ((length(value) == 1L) &&
            is.logical(value) &&
            !is.na(value))
            .CompletionEnv$settings[[what]] <- value
    }
    if (!missing(ops))   checkAndChange(  "ops",   ops)
    if (!missing(ns))    checkAndChange(   "ns",    ns)
    if (!missing(args))  checkAndChange( "args",  args)
    if (!missing(func))  checkAndChange( "func",  func)
    if (!missing(ipck))  checkAndChange( "ipck",  ipck)
    if (!missing(S3))    checkAndChange(   "S3",    S3)
    if (!missing(data))  checkAndChange( "data",  data)
    if (!missing(help))  checkAndChange( "help",  help)
    if (!missing(argdb)) checkAndChange("argdb", argdb)
    if (!missing(files)) checkAndChange("files", files)
    if (!missing(quotes))checkAndChange("quotes", quotes)
    invisible()
}





## modifies options (adapted from similar functions in lattice):

rc.getOption <- function(name)
{
    get("options", envir = .CompletionEnv)[[name]]
}

rc.options <- function(...)
{
    new <- list(...)
    if (is.null(names(new)) && length(new) == 1L && is.list(new[[1L]]))
        new <- new[[1L]]
    old <- .CompletionEnv$options

    ## if no args supplied, returns full options list
    if (length(new) == 0L) return(old)
    ## typically getting options
    nm <- names(new)
    if (is.null(nm)) return(old[unlist(new)])

    isNamed <- nm != ""
    if (any(!isNamed)) nm[!isNamed] <- unlist(new[!isNamed])

    ## so now everything has non-"" names, but only the isNamed ones
    ## should be set.  Everything should be returned though.

    retVal <- old[nm]
    names(retVal) <- nm
    nm <- nm[isNamed]
    .CompletionEnv$options <- modifyList(old, new[nm])
    invisible(retVal)
}



## summarizes results of last completion attempt:

rc.status <- function()
{
    ## eapply(.CompletionEnv, function(x) x, all.names = TRUE)
    as.list(.CompletionEnv)
}


### Everything below is unexported


## accessors called from C (also .completeToken below):

.assignToken         <- function(text)  assign("token",      text,  envir = .CompletionEnv)
.assignLinebuffer    <- function(line)  assign("linebuffer", line,  envir = .CompletionEnv)
.assignStart         <- function(start) assign("start",      start, envir = .CompletionEnv)
.assignEnd           <- function(end)   assign("end",        end,   envir = .CompletionEnv)
.setFileComp         <- function(state) assign("fileName",   state, envir = .CompletionEnv)

.retrieveCompletions <- function()  unique(get("comps",             envir = .CompletionEnv))
.getFileComp         <- function()         get("fileName",          envir = .CompletionEnv)




## The following function is not required for GNU readline, but can be
## useful if one needs to break a line into tokens.  It requires
## linebuffer and end to be already set, and itself sets token and
## start.  It returns the token.

## FIXME: should this use getOption("rl_word_breaks")?

.guessTokenFromLine <-
    function(linebuffer = .CompletionEnv[["linebuffer"]],
             end = .CompletionEnv[["end"]],
             update = TRUE)
    ## update=TRUE changes 'start' and 'token', otherwise they are just returned
{
    linebuffer <- substr(linebuffer, 1L, end) # end is cursor, not necessarily line-end 
    ## special rules apply when we are inside quotes (see fileCompletionPreferred() below)
    insideQuotes <- {
        lbss <- head.default(unlist(strsplit(linebuffer, "")), .CompletionEnv[["end"]])
        ((sum(lbss == "'") %% 2 == 1) ||
         (sum(lbss == '"') %% 2 == 1))
    }
    start <-
        if (insideQuotes)
            ## set 'start' to the location of the last quote
            suppressWarnings(gregexpr("['\"]", linebuffer,
                                      perl = TRUE))[[1L]]
        else
            ##                    things that should not cause breaks
            ##                           _____.^._____
            ##                          /             \
            suppressWarnings(gregexpr("[^\\.\\w:?$@[\\]]+",
                                      linebuffer,
                                      perl = TRUE))[[1L]]
    start <- ## 0-indexed
        if (all(start < 0L)) 0L
        else tail.default(start + attr(start, "match.length"), 1L) - 1L
    token <- substr(linebuffer, start + 1L, end)
    if (update) {
        .CompletionEnv[["start"]] <- start
        .CompletionEnv[["token"]] <- token
        .CompletionEnv[["token"]]
    }
    else list(start = start, token = token)
}






## convert a string to something that escapes special regexp
## characters.  Doesn't have to be perfect, especially for characters
## that would cause breaks or be handled elsewhere.  All we really
## need is to handle ".", so that e.g. "heat." doesn't match
## "heatmap".


makeRegexpSafe <- function(s)
{
    ## the following can cause errors otherwise
    s <- gsub("\\", "\\\\", s, fixed = TRUE) ## has to be the first
    s <- gsub("(", "\\(", s, fixed = TRUE)
    s <- gsub("*", "\\*", s, fixed = TRUE)
    s <- gsub("+", "\\+", s, fixed = TRUE)
    s <- gsub("?", "\\?", s, fixed = TRUE)
    s <- gsub("[", "\\[", s, fixed = TRUE)
    s <- gsub("{", "\\{", s, fixed = TRUE)
    ## s <- gsub("]", "\\]", s, fixed = TRUE) # necessary?
    ## these are wildcards that we want to interpret literally
    s <- gsub(".", "\\.", s, fixed = TRUE)
    s <- gsub("^", "\\^", s, fixed = TRUE)
    ## what else?
    s
}



## Operators that are handled specially.  Order is important, ::: must
## come before :: (because :: will match :::)

specialOps <- c("$", "@", ":::", "::", "?", "[", "[[")


specialOpCompletionsHelper <- function(op, suffix, prefix)
{
    tryToEval <- function(s)
    {
        try(eval(parse(text = s), envir = .GlobalEnv), silent = TRUE)
    }
    switch(op,
           "$" = {
               if (.CompletionEnv$settings[["ops"]])
               {
                   object <- tryToEval(prefix)
                   if (inherits(object, "try-error")) ## nothing else to do
                       suffix
                   else
                   {
                       ## ## suffix must match names(object) (or ls(object) for environments)
                       .DollarNames(object, pattern = sprintf("^%s", makeRegexpSafe(suffix)))
                   }
               } else suffix
           },
           "@" = {
               if (.CompletionEnv$settings[["ops"]])
               {
                   object <- tryToEval(prefix)
                   if (inherits(object, "try-error")) ## nothing else to do
                       suffix
                   else
                   {
                       grep(sprintf("^%s", makeRegexpSafe(suffix)),
                            methods::slotNames(object), value = TRUE)
                   }
               } else suffix
           },
           "::" = {
               if (.CompletionEnv$settings[["ns"]])
               {
                   nse <- try(getNamespaceExports(prefix), silent = TRUE)
                   if (inherits(nse, "try-error")) ## nothing else to do
                       suffix
                   else
                   {
                       grep(sprintf("^%s", makeRegexpSafe(suffix)),
                            nse, value = TRUE)
                   }
               } else suffix
           },
           ":::" = {
               if (.CompletionEnv$settings[["ns"]])
               {
                   ns <- try(getNamespace(prefix), silent = TRUE)
                   if (inherits(ns, "try-error")) ## nothing else to do
                       suffix
                   else
                   {
                       ls(ns,
                          all.names = TRUE,
                          pattern = sprintf("^%s", makeRegexpSafe(suffix)))
                   }
               } else suffix
           },
           "[" = ,  # can't think of anything else to do
           "[[" = {
               comps <- normalCompletions(suffix)
               if (length(comps)) comps
               else suffix
           })
}




specialOpLocs <- function(text)
{
    ## does text contain a special operator?  There may be multiple
    ## occurrences, and we want the last one (whereas regexpr gives
    ## the first). So...

    ge <-
        sapply(specialOps,
               function(s) gregexpr(s, text, fixed = TRUE)[[1L]],
               simplify = FALSE)
    ## this gets the last ones
    ge <- sapply(ge, tail.default, 1)
    ge <- ge[ge > 0]
}



## accessing the help system: should allow anything with an index entry
## this just looks at packages on the search path.

matchAvailableTopics <- function(prefix, text)
{
    .readAliases <- function(path) {
        if(file.exists(f <- file.path(path, "help", "aliases.rds")))
            names(readRDS(f))
        else if(file.exists(f <- file.path(path, "help", "AnIndex")))
            ## aliases.rds was introduced before 2.10.0, as can phase this out
            scan(f, what = list("", ""), sep = "\t", quote = "",
                 na.strings = "", quiet = TRUE)[[1L]]
        else character()
    }
    if (length(text) != 1L || text == "") return (character())
    ## Update list of help topics if necessary
    pkgpaths <- searchpaths()[substr(search(), 1L, 8L) == "package:"]
    if (!identical(basename(pkgpaths), .CompletionEnv[["attached_packages"]])) {
        assign("attached_packages",
               basename(pkgpaths),
               envir = .CompletionEnv)
        assign("help_topics",
               unique(unlist(lapply(pkgpaths, .readAliases))),
               envir = .CompletionEnv)
    }
    aliases <- .CompletionEnv[["help_topics"]]
    ans <- grep(sprintf("^%s", makeRegexpSafe(text)), aliases, value = TRUE)
    if (nzchar(prefix)) {
        tmp <- grep(sprintf("-%s$", prefix), ans, value = TRUE)
        if (length(tmp)) substring(tmp, 1, nchar(tmp) - nchar(prefix) - 1L)
        else character(0)
    }
    else ans
}



## this is for requests of the form ?suffix[TAB] or prefix?suffix[TAB]

## can be improved when prefix is non-trivial, but that is rarely used
## (on the other hand, can be useful for exploring S4 methods; but
## usage of ? needs to be fixed in R first).  Anyway, that case is not
## currently handled

helpCompletions <- function(prefix = "", suffix)
{
    nc <-
        if (.CompletionEnv$settings[["help"]])
            matchAvailableTopics(prefix, suffix)
        else
            normalCompletions(suffix, check.mode = FALSE)
    if (length(nc)) sprintf("%s?%s", prefix, nc)
    else character()
}


specialCompletions <- function(text, spl)
{

    ## we'll only try to complete after the last special operator, and
    ## assume that everything before is meaningfully complete.  A more
    ## sophisticated version of this function may decide to do things
    ## differently.

    ## Note that this will involve evaluations, which may have side
    ## effects.  This (side-effects) would not happen normally (except
    ## of lazy loaded symbols, which most likely would have been
    ## evaluated shortly anyway), because explicit function calls
    ## (with parentheses) are not evaluated.  In any case, these
    ## evaluations won't happen if settings$ops==FALSE

    ## spl (locations of matches) is guaranteed to be non-empty

    wm <- which.max(spl)
    op <- names(spl)[wm]
    opStart <- spl[wm]
    opEnd <- opStart + nchar(op)

    if (opStart < 1) return(character()) # shouldn't happen
    prefix <- substr(text, 1L, opStart - 1L)
    suffix <- substr(text, opEnd, 1000000L)

    if (op == "?") return(helpCompletions(prefix, suffix))

    if (opStart <= 1) return(character()) # not meaningful

    ## ( breaks words, so prefix should not involve function calls,
    ## and thus, hopefully no side-effects.

    comps <- specialOpCompletionsHelper(op, suffix, prefix)
    if (length(comps) == 0L) comps <- ""
    sprintf("%s%s%s", prefix, op, comps)
}



## completions on special keywords (subset of those in gram.c).  Some
## issues with parentheses: e.g. mode(get("repeat")) is "function", so
## it is normally completed with a left-paren appended, but that is
## not normal usage.  Putting it here means that both 'repeat' and
## 'repeat(' will be valid completions (as they should be)


keywordCompletions <- function(text)
{
    grep(sprintf("^%s", makeRegexpSafe(text)),
         c("NULL", "NA", "TRUE", "FALSE", "Inf", "NaN",
           "NA_integer_", "NA_real_", "NA_character_", "NA_complex_",
           "repeat ", "in ", "next ", "break ", "else "),
         value = TRUE)
}




## 'package' environments in the search path.  These will be completed
## with a :: IIRC, that works even if there's no namespace, but I
## haven't actually checked.



attachedPackageCompletions <- function(text, add = rc.getOption("package.suffix"))
{
    if (.CompletionEnv$settings[["ns"]])
    {
        s <- grep("^package", search(), value = TRUE)
        comps <-
            grep(sprintf("^%s", makeRegexpSafe(text)),
                 substr(s, 9L, 1000000L),
                 value = TRUE)
        if (length(comps) && !is.null(add))
            sprintf("%s%s", comps, add)
        else
            comps
    }
    else character()
}



## this provides the most basic completion, looking for completions in
## the search path using apropos, plus keywords.  Plus completion on
## attached packages if settings$ns == TRUE


normalCompletions <-
    function(text, check.mode = TRUE,
             add.fun = rc.getOption("function.suffix"))
{
    ## use apropos or equivalent
    if (text == "") character() ## too many otherwise
    else
    {
        comps <- apropos(sprintf("^%s", makeRegexpSafe(text)), ignore.case = FALSE)
        if (.CompletionEnv$settings[["func"]] && check.mode && !is.null(add.fun))
        {
            which.function <- sapply(comps, function(s) exists(s, mode = "function"))
            if (any(which.function))
                comps[which.function] <-
                    sprintf("%s%s", comps[which.function], add.fun)
            ##sprintf("\033[31m%s\033[0m%s", comps[which.function], add.fun)
        }
        c(comps, keywordCompletions(text), attachedPackageCompletions(text))
    }
}


## completion on function arguments.  This involves the most work (as
## we need to look back in the line buffer to figure out which
## function we are inside, if any), and is also potentially intensive
## when many functions match the function that we are supposedly in
## (which will happen for common generic functions like print (we are
## very optimistic here, erring on the side of
## whatever-the-opposite-of-caution-is (our justification being that
## erring on the side of caution is practically useless and not erring
## at all is expensive to the point of being impossible (we really
## don't want to evaluate the dotplot() call in "print(dotplot(x),
## positi[TAB] )" ))))


## this defines potential function name boundaries

breakRE <- "[^\\.\\w]"
## breakRE <- "[ \t\n \\\" '`><=-%;,&}\\\?\\\+\\\{\\\|\\\(\\\)\\\*]"




## for some special functions like library, data, etc, normal
## completion is rarely meaningful, especially for the first argument.
## Unfortunately, knowing whether the token being completed is the
## first arg of such a function involves more work than we would
## normally want to do.  On the other hand, inFunction() below already
## does most of this work, so we will add a piece of code (mostly
## irrelevant to its primary purpose) to indicate this.  The following
## two functions are just wrappers to access and modify this
## information.


setIsFirstArg <- function(v)
    .CompletionEnv[["isFirstArg"]] <- v

getIsFirstArg <- function()
    .CompletionEnv[["isFirstArg"]]



inFunction <-
    function(line = .CompletionEnv[["linebuffer"]],
             cursor = .CompletionEnv[["start"]])
{
    ## are we inside a function? Yes if the number of ( encountered
    ## going backwards exceeds number of ).  In that case, we would
    ## also like to know what function we are currently inside
    ## (ideally, also what arguments to it have already been supplied,
    ## but let's not dream that far ahead).

    parens <-
        sapply(c("(", ")"),
               function(s) gregexpr(s, substr(line, 1L, cursor), fixed = TRUE)[[1L]],
               simplify = FALSE)
    ## remove -1's
    parens <- lapply(parens, function(x) x[x > 0])

    ## The naive algo is as follows: set counter = 0; go backwards
    ## from cursor, set counter-- when a ) is encountered, and
    ## counter++ when a ( is encountered.  We are inside a function
    ## that starts at the first ( with counter > 0.

    temp <-
        data.frame(i = c(parens[["("]], parens[[")"]]),
                   c = rep(c(1, -1), sapply(parens, length)))
    if (nrow(temp) == 0) return(character())
    temp <- temp[order(-temp$i), , drop = FALSE] ## order backwards
    wp <- which(cumsum(temp$c) > 0)
    if (length(wp)) # inside a function
    {
        ## return guessed name of function, letting someone else
        ## decide what to do with that name

        index <- temp$i[wp[1L]]
        prefix <- substr(line, 1L, index - 1L)
        suffix <- substr(line, index + 1L, cursor + 1L)

        ## note in passing whether we are the first argument (no '='
        ## and no ',' in suffix)

        if ((length(grep("=", suffix, fixed = TRUE)) == 0L) &&
            (length(grep(",", suffix, fixed = TRUE)) == 0L))
            setIsFirstArg(TRUE)


        if ((length(grep("=", suffix, fixed = TRUE))) &&
            (length(grep(",", substr(suffix,
                                     tail.default(gregexpr("=", suffix, fixed = TRUE)[[1L]], 1L),
                                     1000000L), fixed = TRUE)) == 0L))
        {
            ## we are on the wrong side of a = to be an argument, so
            ## we don't care even if we are inside a function
            return(character())
        }
        else ## guess function name
        {
            possible <- suppressWarnings(strsplit(prefix, breakRE, perl = TRUE))[[1L]]
            possible <- possible[possible != ""]
            if (length(possible)) return(tail.default(possible, 1))
            else return(character())
        }
    }
    else # not inside function
    {
        return(character())
    }
}


argNames <-
    function(fname, use.arg.db = .CompletionEnv$settings[["argdb"]])
{
    if (use.arg.db) args <- .FunArgEnv[[fname]]
    if (!is.null(args)) return(args)
    ## else
    args <- do.call(argsAnywhere, list(fname))
    if (is.null(args))
        character()
    else if (is.list(args))
        unlist(lapply(args, function(f) names(formals(f))))
    else
        names(formals(args))
}



specialFunctionArgs <- function(fun, text)
{
    ## certain special functions have special possible arguments.
    ## This is primarily applicable to library and require, for which
    ## rownames(installed.packages()).  This is disabled by default,
    ## because the first call to installed.packages() can be time
    ## consuming, e.g. on a network file system.  However, the results
    ## are cached, so subsequent calls are not that expensive.

    switch(EXPR = fun,

           library = ,
           require = {
               if (.CompletionEnv$settings[["ipck"]])
               {
                   grep(sprintf("^%s", makeRegexpSafe(text)),
                        rownames(installed.packages()), value = TRUE)
               }
               else character()
           },

           data = {
               if (.CompletionEnv$settings[["data"]])
               {
                   grep(sprintf("^%s", makeRegexpSafe(text)),
                        data()$results[, "Item"], value = TRUE)
               }
               else character()
           },

           ## otherwise,
           character())
}



functionArgs <-
    function(fun, text,
             S3methods = .CompletionEnv$settings[["S3"]],
             S4methods = FALSE,
             add.args = rc.getOption("funarg.suffix"))
{
    if (length(fun) < 1L || any(fun == "")) return(character())
    specialFunArgs <- specialFunctionArgs(fun, text)
    if (S3methods && exists(fun, mode = "function"))
        fun <-
            c(fun,
              tryCatch(methods(fun),
                       warning = function(w) {},
                       error = function(e) {}))
    if (S4methods) warning("cannot handle S4 methods yet")
    allArgs <- unique(unlist(lapply(fun, argNames)))
    ans <- grep(sprintf("^%s", makeRegexpSafe(text)), allArgs, value = TRUE)
    if (length(ans) && !is.null(add.args))
        ans <- sprintf("%s%s", ans, add.args)
    c(specialFunArgs, ans)
}



## Note: Inside the C code, we redefine
## rl_attempted_completion_function rather than
## rl_completion_entry_function, which means that if
## retrieveCompletions() returns a length-0 result, by default the
## fallback filename completion mechanism will be used.  This is not
## quite the best way to go, as in certain (most) situations filename
## completion will definitely be inappropriate even if no valid R
## completions are found.  We could return "" as the only completion,
## but that produces an irritating blank line on
## list-possible-completions (or whatever the correct name is).
## Instead (since we don't want to reinvent the wheel), we use the
## following scheme: If the character just preceding our token is " or
## ', we immediately go to file name completion.  If not, we do our
## stuff, and disable file name completion (using
## .Call("RCSuppressFileCompletion")) even if we don't find any
## matches.

## Note that under this scheme, filename completion will fail
## (possibly in unexpected ways) if the partial name contains 'unusual
## characters', namely ones that have been set (see C code) to cause a
## word break because doing so is meaningful in R syntax (e.g. "+",
## "-" ("/" is exempt (and treated specially below) because of its
## ubiquitousness in UNIX file names (where this package is most
## likely to be used))


## decide whether to fall back on filename completion.  Yes if the
## number of quotes between the cursor and the beginning of the line
## is an odd number.

## FIXME: should include backtick (`)? May be useful, but needs more
## thought; e.g., should imply not-filename, but rather variable
## names.  Must cooperate with the if (isInsideQuotes()) branch in
## .completeToken().

isInsideQuotes <-
fileCompletionPreferred <- function()
{
    ((st <- .CompletionEnv[["start"]]) > 0 && {

        ## yes if the number of quote signs to the left is odd
        linebuffer <- .CompletionEnv[["linebuffer"]]
        lbss <- head.default(unlist(strsplit(linebuffer, "")), .CompletionEnv[["end"]])
        ((sum(lbss == "'") %% 2 == 1) ||
         (sum(lbss == '"') %% 2 == 1))

    })
}


## File name completion, used if settings$quotes == TRUE.  Front ends
## that can do filename completion themselves should probably not use
## this if they can do a better job.

correctFilenameToken <- function()
{
    ## Helper function

    ## If a file name contains spaces, the token will only have the
    ## part after the last space.  This function tries to recover the
    ## complete initial part.

    ## Find part between last " or '
    linebuffer <- .CompletionEnv[["linebuffer"]]
    lbss <- head.default(unlist(strsplit(linebuffer, "")), .CompletionEnv[["end"]])
    whichDoubleQuote <- lbss == '"'
    whichSingleQuote <- lbss == "'"
    insideDoubleQuote <- (sum(whichDoubleQuote) %% 2 == 1)
    insideSingleQuote <- (sum(whichSingleQuote) %% 2 == 1)
    loc.start <-
        if (insideDoubleQuote && insideSingleQuote)
        {
            ## Should not happen, but if it does, should take whichever comes later
            max(which(whichDoubleQuote), which(whichSingleQuote))
        }
        else if (insideDoubleQuote)
            max(which(whichDoubleQuote))
        else if (insideSingleQuote)
            max(which(whichSingleQuote))
        else ## should not happen, abort non-intrusively
            .CompletionEnv[["start"]]
    substring(linebuffer, loc.start + 1L, .CompletionEnv[["end"]])
}




fileCompletions <- function(token)
{
    ## uses Sys.glob (conveniently introduced in 2.5.0)

    ## token may not start just after the begin quote, e.g., if spaces
    ## are included.  Get 'correct' partial file name by looking back
    ## to begin quote
    pfilename <- correctFilenameToken()

    ## Sys.glob doesn't work without expansion.  Is that intended?
    pfilename.expanded <- path.expand(pfilename)
    comps <- Sys.glob(sprintf("%s*", pfilename.expanded), dirmark = TRUE)

    ## If there is only one completion (and it's a directory), also
    ## include files inside in list of completions.  This is not
    ## particularly useful, but without this, readline tends to add an
    ## end-quote (if sole completion) which is irritating if one is
    ## actually looking for something inside the directory.  Note that
    ## we don't actually test to see if it's a directory, because if
    ## it is not, list.files() will simply return character(0).
    if (length(comps) == 1 && substring(comps, nchar(comps), nchar(comps)) == "/") {
        filesInside <- list.files(comps, all.files = TRUE, full.names = FALSE, no.. = TRUE)
        if (length(filesInside)) comps <- c(comps, file.path(comps, filesInside))
    }

    ## for things that only extend beyond the cursor, need to
    ## 'unexpand' path
    if (pfilename.expanded != pfilename)
        comps <- sub(path.expand("~"), "~", comps, fixed = TRUE)

    ## for tokens that were non-trivially corrected by adding prefix,
    ## need to delete extra part
    if (pfilename != token)
        comps <- substring(comps, nchar(pfilename) - nchar(token) + 1L, 1000L)
    comps
}





## .completeToken() is the primary interface, and does the actual
## completion when called from C code.


.completeToken <- function()
{
    text <- .CompletionEnv[["token"]]
    if (isInsideQuotes())
    {

        ## If we're in here, that means we think the cursor is inside
        ## quotes.  In most cases, this means that standard filename
        ## completion is more appropriate, but probably not if we're
        ## trying to access things of the form x["foo... or x$"foo...
        ## The following tries to figure this out, but it won't work
        ## in all cases (e.g. x[, "foo<TAB>"])

        ## We assume that whoever determines our token boundaries
        ## considers quote signs as a breaking symbol.


        ## If the 'quotes' setting is FALSE, we will make no attempt to
        ## do filename completion (this is likely to happen with
        ## front-ends that are capable of doing their own file name
        ## completion; such front-ends can fall back to their native
        ## file completion when rc.status("fileName") is TRUE.

        if (.CompletionEnv$settings[["quotes"]])
        {

            ## ## This was used to make a guess whether we are in
            ## ## special situations like ::, ?, [, etc.  But from R
            ## ## 3.0.0 we re-evaluate the token based from the
            ## ## begin-quote, so this is postponed.  This part can be
            ## ## deleted once this is stable enough.
            ## st <- .CompletionEnv[["start"]]
            ## probablyNotFilename <-
            ##     ((st > 2L &&
            ##       ((prequote <- substr(.CompletionEnv[["linebuffer"]], st-1L, st-1L)) %in% c("?", "[", ":", "$"))) ||
            ##      (st == 2L &&
            ##       ((prequote <- substr(.CompletionEnv[["linebuffer"]], st-1L, st-1L)) %in% c("?")))
            ##      )

            ## FIXME|TODO: readline (and maybe other backends) will
            ## usually use a fixed set of breakpoints to detect
            ## tokens.  If we are handling quotes ourselves, the more
            ## likely correct token is everything from the last
            ## unclosed quote onwards (which may include spaces,
            ## punctuations, etc. that would normally cause breaks).
            ## We already do this when we guess the token ourselves
            ## (e.g., for Windows) (and also in the fileCompletions()
            ## call below using correctFilenameToken()), and can
            ## re-use that here.  The problem is that for other
            ## backends a token may already have been determined, and
            ## that's what we will need to use.  We can still fake it
            ## by using the correct token but substracting the extra
            ## part when providing completions, but that will need
            ## some work.

            ## Related to that: if we implement that, should also
            ## check before for '<type>?' and move to help completion
            ## if so.
            
### str(correctFilenameToken())
### str(.guessTokenFromLine(update = FALSE))

            ## TODO: For extra credit, we could also allow for
            ## spaces like in 'package ? grid', but will leave
            ## that for the future (maybe some regexp magic will
            ## make this simple)

            fullToken <- .guessTokenFromLine(update = FALSE)
            probablyHelp <- (fullToken$start >= 2L &&
                             ((substr(.CompletionEnv[["linebuffer"]],
                                      fullToken$start-1L,
                                      fullToken$start-1L)) == "?"))
            if (probablyHelp) {
                fullToken$prefix <- .guessTokenFromLine(end = fullToken$start - 2, update = FALSE)$token
            }
            probablyName <- ((fullToken$start > 2L &&
                              ((substr(.CompletionEnv[["linebuffer"]],
                                       fullToken$start-1L,
                                       fullToken$start-1L)) == "$"))
                             ||
                             (fullToken$start > 3L &&
                              ((substr(.CompletionEnv[["linebuffer"]],
                                       fullToken$start-2L,
                                       fullToken$start-1L)) == "[[")))
            probablyNamespace <- (fullToken$start > 3L &&
                                  ((substr(.CompletionEnv[["linebuffer"]],
                                           fullToken$start-2L,
                                           fullToken$start-1L)) %in% c("::")))
            ## in anticipation that we will handle this eventually:
            probablyBacktick <- (fullToken$start >= 1L &&
                                 ((substr(.CompletionEnv[["linebuffer"]],
                                          fullToken$start,
                                          fullToken$start)) %in% c("`")))
                
            probablySpecial <- probablyHelp || probablyName || probablyNamespace

            ## str(list(probablyHelp = probablyHelp,
            ##          probablyName = probablyName,
            ##          probablyNamespace = probablyNamespace,
            ##          probablyBacktick = probablyBacktick,
            ##          probablySpecial = probablySpecial))

            ## For now, we only handle probablyHelp, and just decline
            ## to do filename completion if any of the other special
            ## situations are detected (but don't try to complete).

            tentativeCompletions <-
                if (probablyHelp) {
                    substring(helpCompletions(fullToken$prefix, fullToken$token),
                              2L + nchar(fullToken$prefix), 1000L)    # drop initial "prefix + ?"
                }
                else if (!probablySpecial)
                    fileCompletions(fullToken$token) # FIXME: but not if probablyBacktick
            .setFileComp(FALSE)
            ## str(c(fullToken, list(comps = tentativeCompletions)))
            ## Adjust for self-computed token
            .CompletionEnv[["comps"]] <-
                substring(tentativeCompletions,
                          1L + nchar(fullToken$token) - nchar(text),
                          1000L)
        }
        else
        {
            .CompletionEnv[["comps"]] <- character()
            .setFileComp(TRUE)
        }
    }
    else
    {
        .setFileComp(FALSE)
        setIsFirstArg(FALSE) # might be changed by inFunction() call
        ## make a guess at what function we are inside
        guessedFunction <-
            if (.CompletionEnv$settings[["args"]])
                inFunction(.CompletionEnv[["linebuffer"]],
                           .CompletionEnv[["start"]])
            else ""
        .CompletionEnv[["fguess"]] <- guessedFunction

        ## if this is not "", then we want to add possible arguments
        ## of that function(s) (methods etc).  Should be character()
        ## if nothing matches

        fargComps <- functionArgs(guessedFunction, text)

        if (getIsFirstArg() && length(guessedFunction) &&
            guessedFunction %in%
            c("library", "require", "data"))
        {
            .CompletionEnv[["comps"]] <- fargComps
            ## don't try anything else
            return()
        }

        ## Is there an arithmetic operator in there in there?  If so,
        ## work on the part after that and append to prefix before
        ## returning.  It would have been easier if these were
        ## word-break characters, but that potentially interferes with
        ## filename completion.

        ## lastArithOp <- tail(gregexpr("/", text, fixed = TRUE)[[1L]], 1)
        lastArithOp <- tail.default(gregexpr("[\"'^/*+-]", text)[[1L]], 1)
        if (haveArithOp <- (lastArithOp > 0))
        {
            prefix <- substr(text, 1L, lastArithOp)
            text <- substr(text, lastArithOp + 1L, 1000000L)
        }

        spl <- specialOpLocs(text)
        comps <-
            if (length(spl))
                specialCompletions(text, spl)
            else
            {
                ## should we append a left-paren for functions?
                ## Usually yes, but not when inside certain special
                ## functions which often take other functions as
                ## arguments

                appendFunctionSuffix <-
                    !any(guessedFunction %in%

                         c("help", "args", "formals", "example",
                           "do.call", "environment", "page", "apply",
                           "sapply", "lapply", "tapply", "mapply",
                           "methods", "fix", "edit"))

                normalCompletions(text, check.mode = appendFunctionSuffix)
            }
        if (haveArithOp && length(comps))
        {
            comps <- paste0(prefix, comps)
        }
        comps <- c(comps, fargComps)
        .CompletionEnv[["comps"]] <- comps
    }
}



## support functions that attempt to provide tools useful specifically
## for the Windows Rgui.


## Note: even though these are unexported functions, changes in the
## API should be noted in man/rcompgen.Rd


.win32consoleCompletion <-
    function(linebuffer, cursorPosition,
             check.repeat = TRUE,
             minlength = -1)
{
    isRepeat <- ## is TAB being pressed repeatedly with this combination?
        if (check.repeat)
            (linebuffer == .CompletionEnv[["linebuffer"]] &&
             cursorPosition == .CompletionEnv[["end"]])
        else TRUE

    .assignLinebuffer(linebuffer)
    .assignEnd(cursorPosition)
    .guessTokenFromLine()
    token <- .CompletionEnv[["token"]]
    comps <-
        if (nchar(token, type = "chars") < minlength) character()
        else
        {
            .completeToken()
            .retrieveCompletions()
        }

    ## FIXME: no idea how much of this is MBCS-safe

    if (length(comps) == 0L)
    {
        ## no completions
        addition <- ""
        possible <- character()
    }
    else if (length(comps) == 1L)
    {
        ## FIXME (maybe): in certain cases the completion may be
        ## shorter than the token (e.g. when trying to complete on an
        ## impossible name inside a list).  It's debatable what the
        ## behaviour should be in this case, but readline and Emacs
        ## actually delete part of the token (at least currently).  To
        ## achieve this in Rgui one would need to do somewhat more
        ## work than I'm ready to do right now (especially since it's
        ## not clear that this is the right thing to do to begin
        ## with).  So, in this case, I'll just pretend that no
        ## completion was found.

        addition <- substr(comps, nchar(token, type = "chars") + 1L, 100000L)
        possible <- character()
    }
    else if (length(comps) > 1L)
    {
        ## more than one completion.  The right thing to is to extend
        ## the line by the unique part if any, and list the multiple
        ## possibilities otherwise.

        additions <- substr(comps, nchar(token, type = "chars") + 1L, 100000L)
        if (length(table(substr(additions, 1L, 1L))) > 1L)
        {
            ## no unique substring
            addition <- ""
            possible <-
                if (isRepeat) capture.output(cat(format(comps, justify = "left"), fill = TRUE))
                else character()
        }
        else
        {
            ## need to figure out maximal unique substr
            keepUpto <- 1
            while (length(table(substr(additions, 1L, keepUpto))) == 1L)
                keepUpto <- keepUpto + 1L
            addition <- substr(additions[1L], 1L, keepUpto - 1L)
            possible <- character()
        }
    }
    list(addition = addition,
         possible = possible,
         comps = paste(comps, collapse = " "))
}




## usage:

## .addFunctionInfo(foo = c("arg1", "arg2"), bar = c("a", "b"))

.addFunctionInfo <- function(...)
{
    dots <- list(...)
    for (nm in names(dots))
        .FunArgEnv[[nm]] <- dots[[nm]]
}

.initialize.argdb <-
    function()
{
    ## lattice

    lattice.common <-
        c("data", "allow.multiple", "outer", "auto.key", "aspect",
          "panel", "prepanel", "scales", "strip", "groups", "xlab",
          "xlim", "ylab", "ylim", "drop.unused.levels", "...",
          "default.scales", "subscripts", "subset", "formula", "cond",
          "aspect", "as.table", "between", "key", "legend", "page",
          "main", "sub", "par.strip.text", "layout", "skip", "strip",
          "strip.left", "xlab.default", "ylab.default", "xlab",
          "ylab", "panel", "xscale.components", "yscale.components",
          "axis", "index.cond", "perm.cond", "...", "par.settings",
          "plot.args", "lattice.options")

    densityplot <-
        c("plot.points", "ref", "groups", "jitter.amount",
          "bw", "adjust", "kernel", "weights", "window", "width",
          "give.Rkern", "n", "from", "to", "cut", "na.rm")

    panel.xyplot <-
        c("type", "groups", "pch", "col", "col.line",
          "col.symbol", "font", "fontfamily", "fontface", "lty",
          "cex", "fill", "lwd", "horizontal")

    .addFunctionInfo(xyplot.formula = c(lattice.common, panel.xyplot),
                     densityplot.formula = c(lattice.common, densityplot))

    ## grid

    grid.clip <-
        c("x", "y", "width", "height", "just", "hjust", "vjust",
          "default.units", "name", "vp")
    grid.curve <-
        c("x1", "y1", "x2", "y2", "default.units", "curvature",
          "angle", "ncp", "shape", "square", "squareShape", "inflect",
          "arrow", "open", "debug", "name", "gp", "vp")
    grid.polyline <-
        c("x", "y", "id", "id.lengths", "default.units", "arrow",
          "name", "gp", "vp")
    grid.xspline <-
        c("x", "y", "id", "id.lengths", "default.units", "shape",
          "open", "arrow", "repEnds", "name", "gp", "vp")

    .addFunctionInfo(grid.clip = grid.clip,
                     grid.curve = grid.curve,
                     grid.polyline = grid.polyline,
                     grid.xspline = grid.xspline)

    ## par, options

    par <-
        c("xlog", "ylog", "adj", "ann", "ask", "bg", "bty", "cex",
          "cex.axis", "cex.lab", "cex.main", "cex.sub", "cin", "col",
          "col.axis", "col.lab", "col.main", "col.sub", "cra", "crt",
          "csi", "cxy", "din", "err", "family", "fg", "fig", "fin",
          "font", "font.axis", "font.lab", "font.main", "font.sub",
          "gamma", "lab", "las", "lend", "lheight", "ljoin", "lmitre",
          "lty", "lwd", "mai", "mar", "mex", "mfcol", "mfg", "mfrow",
          "mgp", "mkh", "new", "oma", "omd", "omi", "pch", "pin",
          "plt", "ps", "pty", "smo", "srt", "tck", "tcl", "usr",
          "xaxp", "xaxs", "xaxt", "xpd", "yaxp", "yaxs", "yaxt")

    options <- c("add.smooth", "browser", "check.bounds", "continue",
	"contrasts", "defaultPackages", "demo.ask", "device",
	"digits", "dvipscmd", "echo", "editor", "encoding",
	"example.ask", "expressions", "help.search.types",
	"help.try.all.packages", "htmlhelp", "HTTPUserAgent",
	"internet.info", "keep.source", "keep.source.pkgs",
	"locatorBell", "mailer", "max.print", "menu.graphics",
	"na.action", "OutDec", "pager", "papersize",
	"par.ask.default", "pdfviewer", "pkgType", "printcmd",
	"prompt", "repos", "scipen", "show.coef.Pvalues",
	"show.error.messages", "show.signif.stars", "str",
	"stringsAsFactors", "timeout", "ts.eps", "ts.S.compat",
	"unzip", "verbose", "warn", "warning.length", "width")

    .addFunctionInfo(par = par, options = options)

    ## read.csv etc (... passed to read.table)

}



.CompletionEnv <- new.env(hash = FALSE)

## needed to save some overhead in .win32consoleCompletion
assign("linebuffer", "", env = .CompletionEnv)
assign("end", 1, env = .CompletionEnv)

assign("settings",
       list(ops = TRUE, ns = TRUE,
            args = TRUE, func = FALSE,
            ipck = FALSE, S3 = TRUE, data = TRUE,
            help = TRUE, argdb = TRUE,
            files = TRUE, # FIXME: deprecate in favour of quotes
            quotes = TRUE),
       env = .CompletionEnv)

assign("options",
       list(package.suffix = "::",
            funarg.suffix = "=",
            function.suffix = "("),
       env = .CompletionEnv)

## These keeps track of attached packages and available help topics.
## Needs updating only when packages are attached.
assign("attached_packages", character(0), env = .CompletionEnv)
assign("help_topics", character(0), env = .CompletionEnv)


.FunArgEnv <- new.env(hash = TRUE, parent = emptyenv())

.initialize.argdb()

#  File src/library/utils/R/data.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

data <-
function(..., list = character(), package = NULL, lib.loc = NULL,
         verbose = getOption("verbose"), envir = .GlobalEnv)
{
    fileExt <- function(x) {
        db <- grepl("\\.[^.]+\\.(gz|bz2|xz)$", x)
        ans <- sub(".*\\.", "", x)
        ans[db] <-  sub(".*\\.([^.]+\\.)(gz|bz2|xz)$", "\\1\\2", x[db])
        ans
    }

    names <- c(as.character(substitute(list(...))[-1L]), list)

    ## Find the directories of the given packages and maybe the working
    ## directory.
    if(!is.null(package)) {
        if(!is.character(package))
            stop("'package' must be a character string or NULL")
        if(any(package %in% "base"))
            warning("datasets have been moved from package 'base' to package 'datasets'")
        if(any(package %in% "stats"))
           warning("datasets have been moved from package 'stats' to package 'datasets'")
        package[package %in% c("base", "stats")] <- "datasets"
    }
    paths <- find.package(package, lib.loc, verbose = verbose)
    if(is.null(lib.loc))
        paths <- c(path.package(package, TRUE),
                   if(!length(package)) getwd(), # ignored if NULL
                   paths)
    paths <- unique(paths[file.exists(paths)])

    ## Find the directories with a 'data' subdirectory.
    paths <- paths[file_test("-d", file.path(paths, "data"))]

    dataExts <- tools:::.make_file_exts("data")

    if(length(names) == 0L) {
        ## List all possible data sets.

        ## Build the data db.
        db <- matrix(character(), nrow = 0L, ncol = 4L)
        for(path in paths) {
            entries <- NULL
            ## Use "." as the 'package name' of the working directory.
            packageName <-
                if(file_test("-f", file.path(path, "DESCRIPTION")))
                    basename(path)
                else
                    "."
            ## Check for new-style 'Meta/data.rds'
            if(file_test("-f", INDEX <- file.path(path, "Meta", "data.rds"))) {
                entries <- readRDS(INDEX)
            } else {
                ## No index: should only be true for ./data >= 2.0.0
                dataDir <- file.path(path, "data")
                entries <- tools::list_files_with_type(dataDir, "data")
                if(length(entries)) {
                    entries <-
                        unique(tools::file_path_sans_ext(basename(entries)))
                    entries <- cbind(entries, "")
                }
            }
            if(NROW(entries)) {
                if(is.matrix(entries) && ncol(entries) == 2L)
                    db <- rbind(db, cbind(packageName, dirname(path), entries))
                else
                    warning(gettextf("data index for package %s is invalid and will be ignored",
                                     sQuote(packageName)),
                            domain=NA, call.=FALSE)
            }
        }
        colnames(db) <- c("Package", "LibPath", "Item", "Title")

        footer <- if(missing(package))
            paste0("Use ",
                   sQuote(paste("data(package =",
                                ".packages(all.available = TRUE))")),
                   "\n",
                   "to list the data sets in all *available* packages.")
        else
            NULL
        y <- list(title = "Data sets", header = NULL, results = db,
                  footer = footer)
        class(y) <- "packageIQR"
        return(y)
    }

    paths <- file.path(paths, "data")
    for(name in names) {
        found <- FALSE
        for(p in paths) {
            ## does this package have "Rdata" databases?
            if(file_test("-f", file.path(p, "Rdata.rds"))) {
                rds <- readRDS(file.path(p, "Rdata.rds"))
                if(name %in% names(rds)) {
                    ## found it, so copy objects from database
                    found <- TRUE
                    if(verbose)
                        message(sprintf("name=%s:\t found in Rdata.rds", name),
                                domain=NA)
                    thispkg <- sub(".*/([^/]*)/data$", "\\1", p)
                    thispkg <- sub("_.*$", "", thispkg) # versioned installs.
                    thispkg <- paste0("package:", thispkg)
                    objs <- rds[[name]] # guaranteed an exact match
                    lazyLoad(file.path(p, "Rdata"), envir = envir,
                             filter = function(x) x %in% objs)
                    break
		} else if(verbose)
		    message(sprintf("name=%s:\t NOT found in names() of Rdata.rds, i.e.,\n\t%s\n",
				    name, paste(names(rds), collapse=",")),
				domain=NA)
            }
            ## check for zipped data dir
            if(file_test("-f", file.path(p, "Rdata.zip"))) {
                warning("zipped data found for package ",
                        sQuote(basename(dirname(p))),
                        ".\nThat is defunct, so please re-install the package.",
                        domain = NA)
                if(file_test("-f", fp <- file.path(p, "filelist")))
                    files <- file.path(p, scan(fp, what="", quiet = TRUE))
                else {
                    warning(gettextf("file 'filelist' is missing for directory %s", sQuote(p)), domain = NA)
                    next
                }
            } else {
                files <- list.files(p, full.names = TRUE)
            }
            files <- files[grep(name, files, fixed = TRUE)]
            if(length(files) > 1L) {
                ## more than one candidate
                o <- match(fileExt(files), dataExts, nomatch = 100L)
                paths0 <- dirname(files)
		## Next line seems unnecessary to MM (FIXME?)
		paths0 <- factor(paths0, levels = unique(paths0))
                files <- files[order(paths0, o)]
            }
            if(length(files)) {
                ## have a plausible candidate (or more)
                for(file in files) {
                    if(verbose)
                        message("name=", name, ":\t file= ...",
                                .Platform$file.sep, basename(file), "::\t",
                                appendLF = FALSE, domain = NA)
                    ext <- fileExt(file)
                    ## make sure the match is really for 'name.ext'
                    if(basename(file) != paste0(name, ".", ext))
                        found <- FALSE
                    else {
                        found <- TRUE
                        zfile <- file
                        zipname <- file.path(dirname(file), "Rdata.zip")
                        if(file.exists(zipname)) {
                            Rdatadir <- tempfile("Rdata")
                            dir.create(Rdatadir, showWarnings=FALSE)
                            topic <- basename(file)
                            rc <- .External(C_unzip, zipname, topic, Rdatadir, FALSE, TRUE, FALSE, FALSE)
                            if(rc == 0L) zfile <- file.path(Rdatadir, topic)
                        }
                        if(zfile != file) on.exit(unlink(zfile))
                        switch(ext,
                               R = , r = {
                                   ## ensure utils is visible
                                   library("utils")
                                   sys.source(zfile, chdir = TRUE,
                                              envir = envir)
                               },
                               RData = , rdata = , rda =
                               load(zfile, envir = envir),
                               TXT = , txt = , tab = ,
                               tab.gz = , tab.bz2 = , tab.xz = ,
                               txt.gz = , txt.bz2 = , txt.xz =
                               assign(name,
                                      ## ensure default for as.is has not been
                                      ## overridden by options(stringsAsFactor)
                                      read.table(zfile, header = TRUE, as.is = FALSE),
                                      envir = envir),
                               CSV = , csv = ,
                               csv.gz = , csv.bz2 = , csv.xz =
                               assign(name,
                                      read.table(zfile, header = TRUE,
                                                 sep = ";", as.is = FALSE),
                                      envir = envir),
                               found <- FALSE)
                    }
                    if (found) break # from files
                }
                if(verbose) message(if(!found) "*NOT* ", "found", domain = NA)
            }
            if (found) break # from paths
        }

        if(!found)
            warning(gettextf("data set %s not found", sQuote(name)),
                    domain = NA)
    }
    invisible(names)
}
#  File src/library/utils/R/databrowser.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

browseEnv <- function(envir = .GlobalEnv, pattern,
                      excludepatt = "^last\\.warning",
		      html = .Platform$GUI != "AQUA",
		      expanded = TRUE, properties = NULL,
		      main = NULL, debugMe = FALSE)
{
    objlist <- ls(envir = envir, pattern = pattern)#, all.names = FALSE
    if(length(iX <- grep(excludepatt, objlist)))
        objlist <- objlist[ - iX]
    if(debugMe) { cat("envir= "); print(envir)
		  cat("objlist =\n"); print(objlist) }
    n <- length(objlist)
    if(n == 0L) {
	cat("Empty environment, nothing to do!\n")
	return(invisible())
    }

    str1 <- function(obj) {
	md <- mode(obj)
	lg <- length(obj)
	objdim <- dim(obj)
	if(length(objdim) == 0L)
	    dim.field <- paste("length:", lg)
	else{
	    dim.field <- "dim:"
	    for(i in seq_along(objdim))
		dim.field <- paste(dim.field,objdim[i])
	    if(is.matrix(obj))
		md <- "matrix"
	}
	obj.class <- oldClass(obj)
	if(!is.null(obj.class)) {
	    md <- obj.class[1L]
	    if(inherits(obj, "factor"))
		dim.field <- paste("levels:",length(levels(obj)))
	}
	list(type = md, dim.field = dim.field)
    }

    N <- 0L
    M <- n
    IDS <- rep.int(NA,n)
    NAMES <- rep.int(NA,n)
    TYPES <- rep.int(NA,n)
    DIMS <- rep.int(NA,n)

    IsRoot <- rep.int(TRUE,n)
    Container <- rep.int(FALSE,n)
    ItemsPerContainer <- rep.int(0,n)
    ParentID <- rep.int(-1,n)

    for( objNam in objlist ){
	N <- N+1L
	if(debugMe) cat("  ", N,":", objNam)
	obj    <- get(objNam, envir = envir)

	sOb <- str1(obj)

	if(debugMe) cat(", type=", sOb$type,",", sOb$dim.field,"\n")

	## Fixme : put these 4 in a matrix or data.frame row:
	IDS[N] <- N
	NAMES[N] <- objNam
	TYPES[N] <- sOb$type
	DIMS[N] <- sOb$dim.field

	if(is.recursive(obj) && !is.function(obj) && !is.environment(obj)
	    ## includes "list", "expression", also "data.frame", ..
	   && (lg <- length(obj))) {
	    Container[N] <- TRUE
	    ItemsPerContainer[N] <- lg
	    nm <- names(obj)
	    if(is.null(nm)) nm <- paste0("[[", format(1L:lg), "]]")
	    for(i in 1L:lg) {
		M <- M+1
		ParentID[M] <- N
		if(nm[i] == "") nm[i] <- paste0("[[", i, "]]")

		s.l <- str1(obj[[i]])
		##cat("	   objname:",nm[i],", type=",md.l,",",dim.field.l,"\n")
		IDS   <- c(IDS,M)
		NAMES <- c(NAMES, nm[i])
		TYPES <- c(TYPES, s.l$type)
		DIMS  <- c(DIMS,  s.l$dim.field)
	    }
	}## recursive

	else if(!is.null(class(obj))) {
	    ## treat some special __non-recursive__ classes:
	    if(inherits(obj, "table")) {
		obj.nms <- attr(obj,"dimnames")
		lg <- length(obj.nms)
		if(length(names(obj.nms)) >0)
		    nm <- names(obj.nms)
		else
		    nm <- rep.int("", lg)
		Container[N] <- TRUE
		ItemsPerContainer[N] <- lg
		for(i in 1L:lg){
		    M <- M+1L
		    ParentID[M] <- N
		    if(nm[i] == "") nm[i] = paste0("[[",i,"]]")
		    md.l  <- mode(obj.nms[[i]])
		    objdim.l <- dim(obj.nms[[i]])
		    if(length(objdim.l) == 0L)
			dim.field.l <- paste("length:", length(obj.nms[[i]]))
		    else{
			dim.field.l <- "dim:"
			for(j in seq_along(objdim.l))
			    dim.field.l <- paste(dim.field.l,objdim.l[i])
		    }
		    ##cat("    objname:",nm[i],", type=",md.l,",",dim.field.l,"\n")
		    IDS <- c(IDS,M)
		    NAMES <- c(NAMES, nm[i])
		    TYPES <- c(TYPES, md.l)
		    DIMS <- c(DIMS,dim.field.l)
		}
	    }## "table"

	    else if(inherits(obj, "mts")) {

		nm <- dimnames(obj)[[2L]]
		lg <- length(nm)
		Container[N] <- TRUE
		ItemsPerContainer[N] <- lg
		for(i in 1L:lg){
		    M <- M+1L
		    ParentID[M] <- N
		    md.l  <- mode(obj[[i]])
		    dim.field.l <- paste("length:",dim(obj)[1L])
		    md.l <- "ts"
		    ##cat("    tseries:",nm[i],", type=",md.l,",",dim.field.l,"\n")
		    IDS <- c(IDS,M)
		    NAMES <- c(NAMES, nm[i])
		    TYPES <- c(TYPES, md.l)
		    DIMS <- c(DIMS,dim.field.l)
		}
	    }## "mts"

	} ## recursive or classed

    } ## "for each object"

    if(debugMe) cat(" __end {for}\n ")##; browser()

    Container	      <- c(Container,	  rep.int(FALSE, M-N))
    IsRoot	      <- c(IsRoot,	  rep.int(FALSE, M-N))
    ItemsPerContainer <- c(ItemsPerContainer, rep.int(0, M-N))

    if(is.null(main))
	main <- paste("R objects in", deparse(substitute(envir)))
    if(is.null(properties)) {
	properties <- as.list(c(date = format(Sys.time(), "%Y-%b-%d %H:%M"),
				local({
				    si <- Sys.info()
				    si[c("user","nodename","sysname")]})))
    }
    if(html)
	wsbrowser(IDS, IsRoot, Container, ItemsPerContainer, ParentID,
		  NAMES, TYPES, DIMS, kind = "HTML", main = main,
                  properties = properties, expanded)
    else if(.Platform$GUI == "AQUA") {
        awsbrowser <- get("wsbrowser", envir = as.environment("tools:RGUI"))
 	awsbrowser(as.integer(IDS), IsRoot, Container,
                   as.integer(ItemsPerContainer), as.integer(ParentID),
                   NAMES, TYPES, DIMS)
   } else stop("only 'html = TRUE' is supported on this platform")
}

wsbrowser <- function(IDS, IsRoot, IsContainer, ItemsPerContainer,
		      ParentID, NAMES, TYPES, DIMS, expanded=TRUE,
		      kind = "HTML",
		      main = "R Workspace", properties = list(),
		      browser = getOption("browser"))
{
    if(kind != "HTML")
        stop(gettextf("kind '%s' not yet implemented", kind), domain = NA)

    bold <- function(ch) paste0("<b>",ch,"</b>")
    ital <- function(ch) paste0("<i>",ch,"</i>")
    entry <- function(ch) paste0("<td>",ch,"</td>")
    Par	<- function(ch) paste0("<P>",ch,"</P>")
    Trow <- function(N, ...) {
	if(length(list(...)) != N) stop("wrong number of table row entries")
	paste("<tr>", ..., "</tr>\n")
    }
    catRow <- function(...) cat(Trow(nCol, ...), file = Hfile)

#    n <- length(IDS)
    RootItems <- which(IsRoot)
    NumOfRoots <- length(RootItems)

    props <- properties
    if(length(props)) { ## translate named list into 2-column (vertical) table
	nms <- names(props)
	nms <- unlist(lapply(unlist(lapply(paste0(nms,":"),
					   bold)),
			     entry))
	props <- unlist(lapply(props, entry))
	props <-
	    paste("<table border=2>",
		  paste(Trow(1, paste(nms, props)), collapse=""),
		  "</table>", sep = "\n")
    }
    fname <- file.path(tempdir(), "wsbrowser.html")
    Hfile <- file(fname,"w")

    cat("<html>\n<title>", main, "browser</title>\n<body>",
	"<H1>",main,"</H1>\n",
	if(is.character(props)) Par(props),
	"<table border=1>\n", file = Hfile)
    nCol <- if(expanded) 4L else 3L
    catRow(entry(bold("Object")),
	   if(expanded) entry(bold(ital("(components)"))),
	   entry(bold("Type")),
	   entry(bold("Property")))

    for(i in 1L:NumOfRoots) {
	iid <- RootItems[i]
	catRow(entry(NAMES[iid]),
	       if(expanded) entry(""),
	       entry(ital(TYPES[iid])),
	       entry(DIMS[iid]))
	if(IsContainer[i] && expanded) {
	    items <- which(ParentID == i)
	    for(j in 1L:ItemsPerContainer[i]) {
		id <- IDS[items[j]]
		catRow(entry(""),
		       entry(NAMES[id]),#was paste0("$",NAMES[id]) : ugly for [[i]]
		       entry(ital(TYPES[id])),
		       entry(DIMS[id]))
	    }
	}
    }
    cat("</table>\n</body></html>",file=Hfile)
    close(Hfile)

    switch(.Platform$OS.type,
	   windows = , ## do we need anything here?
	   unix = { url <- fname },
	   )
    if(substr(url, 1L, 1L) != "/")
	url <- paste0("/", url)
    url <- paste0("file://", URLencode(url))

    browseURL(url = url, browser = browser)
    cat(main, "environment is shown in browser",
	if(is.character(browser)) sQuote(browser),"\n")

    invisible(fname)
}
#  File src/library/utils/R/de.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

de.ncols <- function(inlist)
{
    ncols <- matrix(0, nrow=length(inlist), ncol=2L)
    i <- 1L
    for( telt in inlist ) {
	if( is.matrix(telt) ) {
	    ncols[i, 1L] <- ncol(telt)
	    ncols[i, 2L] <- 2L
	}
	else if( is.list(telt) ) {
	    for( telt2 in telt )
		if( !is.vector(telt2) ) stop("wrong argument to 'dataentry'")
	    ncols[i, 1L] <- length(telt)
	    ncols[i, 2L] <- 3L
	}
	else if( is.vector(telt) ) {
	    ncols[i, 1L] <- 1L
	    ncols[i, 2L] <- 1L
	}
	else stop("wrong argument to 'dataentry'")
	i <- i+1L
    }
    return(ncols)
}

de.setup <- function(ilist, list.names, incols)
{
    ilen <- sum(incols)
    ivec <- vector("list", ilen)
    inames <- vector("list", ilen)
    i <- 1L
    k <- 0L
    for( telt in ilist ) {
	k <- k+1L
	if( is.list(telt) ) {
	    y <- names(telt)
	    for( j in seq_along(telt) ) {
		ivec[[i]] <- telt[[j]]
		if( is.null(y) || y[j]=="" )
		    inames[[i]] <- paste0("var", i)
		else inames[[i]] <- y[j]
		i <- i+1L
	    }
	}
	else if( is.vector(telt) ) {
	    ivec[[i]] <- telt
	    inames[[i]] <- list.names[[k]]
	    i <- i+1
	}
	else if( is.matrix(telt) ) {
	    y <- dimnames(telt)[[2L]]
	    for( j in 1L:ncol(telt) ) {
		ivec[[i]] <- telt[, j]
		if( is.null(y) || y[j]=="" )
		    inames[[i]] <- paste0("var", i)
		else inames[[i]] <- y[j]
		i <- i+1L
	    }
	}
	else stop("wrong argument to 'dataentry'")
    }
    names(ivec) <- inames
    return(ivec)
}

de.restore <- function(inlist, ncols, coltypes, argnames, args)
{
    ## take the data in inlist and restore it
    ## to the format described by ncols and coltypes
    p <- length(ncols)
    rlist <- vector("list", length=p)
    rnames <- vector("character", length=p)
    j <- 1L
    lnames <- names(inlist)
    if(p) for(i in 1L:p) {
	if(coltypes[i]==2) {
	    tlen <- length(inlist[[j]])
	    x <- matrix(0, nrow=tlen, ncol=ncols[i])
	    cnames <- vector("character", ncol(x))
	    for( ind1 in 1L:ncols[i]) {
		if(tlen != length(inlist[[j]]) ) {
		    warning("could not restore type information")
		    return(inlist)
		}
		x[, ind1] <- inlist[[j]]
		cnames[ind1] <- lnames[j]
		j <- j+1L
	    }
	    if( nrow(x) == nrow(args[[i]]) )
		rn <- dimnames(args[[i]])[[1L]]
	    else rn <- NULL
	    if( any(cnames!="") )
		dimnames(x) <- list(rn, cnames)
	    rlist[[i]] <- x
	    rnames[i] <- argnames[i]
	}
	else if(coltypes[i]==3) {
	    x <- vector("list", length=ncols[i])
	    cnames <- vector("character", ncols[i])
	    for( ind1 in 1L:ncols[i]) {
		x[[ind1]] <- inlist[[j]]
		cnames[ind1] <- lnames[j]
		j <- j+1L
	    }
	    if( any(cnames!="") )
		names(x) <- cnames
	    rlist[[i]] <- x
	    rnames[i] <- argnames[i]
	}
	else {
	    rlist[[i]] <- inlist[[j]]
	    j <- j+1
	    rnames[i] <- argnames[i]
	}
    }
    names(rlist) <- rnames
    return(rlist)
}

de <- function(..., Modes=list(), Names=NULL)
{
    sdata <- list(...)
    snames <- as.character(substitute(list(...))[-1L])
    if( is.null(sdata) ) {
	if( is.null(Names) ) {
	    odata <- vector("list", length=max(1,length(Modes)))
	}
	else {
	    if( (length(Names) != length(Modes)) && length(Modes) ) {
		warning("'modes' argument ignored")
		Modes <- list()
	    }
	    odata <- vector("list", length=length(Names))
	    names(odata) <- Names
	}
	ncols <- rep.int(1, length(odata))
	coltypes <- rep.int(1, length(odata))
    }
    else {
	ncols <- de.ncols(sdata)
	coltypes <- ncols[, 2L]
	ncols <- ncols[, 1]
	odata <- de.setup(sdata, snames, ncols)
	if(length(Names))
	    if( length(Names) != length(odata) )
		warning("'names' argument ignored")
	    else names(odata) <- Names
	if(length(Modes))
	    if(length(Modes) != length(odata)) {
		warning("'modes' argument ignored")
		Modes <- list()
	    }
    }
    rdata <- dataentry(odata, as.list(Modes))

    if(any(coltypes != 1L)) {
	if(length(rdata) == sum(ncols))
	    rdata <- de.restore(rdata, ncols, coltypes, snames, sdata)
	else warning("could not restore variables properly")
    }
    return(rdata)
}

data.entry <- function(..., Modes=NULL, Names=NULL)
{
    tmp1 <- de(..., Modes=Modes, Names=Names)
    j <- 1L
    nn <- names(tmp1)
    for(i in nn) {
	assign(i, tmp1[[j]], envir=.GlobalEnv)
	j <- j+1L
    }
    if(j == 1L) warning("did not assign() anything")
    invisible(nn)
}
#  File src/library/utils/R/debugger.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

dump.frames <- function(dumpto = "last.dump", to.file = FALSE)
{
    calls <- sys.calls()
    last.dump <- sys.frames()
    names(last.dump) <- limitedLabels(calls)
    last.dump <- last.dump[-length(last.dump)] # remove this function
    attr(last.dump, "error.message") <- geterrmessage()
    class(last.dump) <- "dump.frames"
    if(dumpto != "last.dump") assign(dumpto, last.dump)
    if (to.file) # compress=TRUE is now the default.
        save(list=dumpto, file = paste(dumpto, "rda", sep = "."))
    else assign(dumpto, last.dump, envir=.GlobalEnv)
    invisible()
}

debugger <- function(dump = last.dump)
{
    debugger.look <- function(.selection)
    {
        ## allow e.g. '...' to fail
        for(.obj in ls(envir=dump[[.selection]], all.names=TRUE))
            tryCatch(assign(.obj, get(.obj, envir=dump[[.selection]])),
                     error=function(e) {})
        cat(gettext("Browsing in the environment with call:\n   "),
            calls[.selection], "\n", sep = "")
        rm(.obj, .selection)
        browser()
    }
    if (!inherits(dump, "dump.frames")) {
        cat(gettextf("'dump' is not an object of class %s\n",
                     dQuote("dump.frames")))
        return(invisible())
    }
    err.action <- getOption("error")
    on.exit(options(error=err.action))
    if (length(msg <- attr(dump, "error.message")))
        cat(gettext("Message: "), msg)
    n <- length(dump)
    calls <- names(dump)
    repeat {
        cat(gettext("Available environments had calls:\n"))
        cat(paste0(1L:n, ": ", calls), sep = "\n")
        cat(gettext("\nEnter an environment number, or 0 to exit  "))
        repeat {
            ind <- .Call(C_menu, as.character(calls))
            if(ind <= n) break
        }
        if(ind == 0L) return(invisible())
        debugger.look(ind)
    }
}

## allow for the numbering by menu here
limitedLabels <- function(value, maxwidth = getOption("width") - 5L)
{
    srcrefs <- sapply(value, function(v)
                      if (!is.null(srcref <- attr(v, "srcref"))) {
                          srcfile <- attr(srcref, "srcfile")
                          paste0(basename(srcfile$filename), "#", srcref[1L],": ")
                      } else "")
    value <- paste0(srcrefs, as.character(value))
    if(is.null(maxwidth) || maxwidth < 40L) maxwidth <- 40L
    maxwidth <- min(maxwidth, 1000L)
    strtrim(value, maxwidth)
}

recover <-
  function()
{
    if(.isMethodsDispatchOn()) {
        ## turn off tracing
        tState <- tracingState(FALSE)
        on.exit(tracingState(tState))
    }
    ## find an interesting environment to start from
    calls <- sys.calls()
    from <- 0L
    n <- length(calls)
    if(identical(sys.function(n), recover))
        ## options(error=recover) produces a call to this function as an object
        n <- n - 1L
    ## look for a call inserted by trace() (and don't show frames below)
    ## this level.
    for(i in rev(seq_len(n))) {
        calli <- calls[[i]]
        fname <- calli[[1L]]
        ## deparse can use more than one line
        if(!is.na(match(deparse(fname)[1L],
                        c("methods::.doTrace", ".doTrace")))) {
            from <- i-1L
            break
        }
    }
  ## if no trace, look for the first frame from the bottom that is not
    ## stop or recover
    if(from == 0L)
      for(i in rev(seq_len(n))) {
        calli <- calls[[i]]
        fname <- calli[[1L]]
        if(!is.name(fname) ||
           is.na(match(as.character(fname), c("recover", "stop", "Stop")))) {
            from <- i
            break
        }
    }
    if(from > 0L) {
        if(!interactive()) {
            try(dump.frames())
            cat(gettext("recover called non-interactively; frames dumped, use debugger() to view\n"))
            return(NULL)
        }
        else if(identical(getOption("show.error.messages"), FALSE)) # from try(silent=TRUE)?
            return(NULL)
        calls <- limitedLabels(calls[1L:from])
        repeat {
            which <- menu(calls,
                          title="\nEnter a frame number, or 0 to exit  ")
            if(which)
                eval(substitute(browser(skipCalls=skip),
                                list(skip=7-which)), envir = sys.frame(which))
            else
                break
        }
    }
    else
        cat(gettext("No suitable frames for recover()\n"))
}
#  File src/library/utils/R/demo.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

demo <-
function(topic, package = NULL, lib.loc = NULL,
	 character.only = FALSE, verbose = getOption("verbose"),
	 echo = TRUE, ask = getOption("demo.ask"),
         encoding = getOption("encoding"))
{
    paths <- find.package(package, lib.loc, verbose = verbose)

    ## Find the directories with a 'demo' subdirectory.
    paths <- paths[file_test("-d", file.path(paths, "demo"))]
    ## Earlier versions remembered given packages with no 'demo'
    ## subdirectory, and warned about them.

    if(missing(topic)) {
	## List all possible demos.

	## Build the demo db.
	db <- matrix(character(), nrow = 0L, ncol = 4L)
	for(path in paths) {
	    entries <- NULL
	    ## Check for new-style 'Meta/demo.rds', then for '00Index'.
	    if(file_test("-f", INDEX <- file.path(path, "Meta", "demo.rds"))) {
		entries <- readRDS(INDEX)
	    }
	    if(NROW(entries)) {
		db <- rbind(db,
			    cbind(basename(path), dirname(path),
				  entries))
	    }
	}
	colnames(db) <- c("Package", "LibPath", "Item", "Title")

	footer <- if(missing(package))
	    paste0("Use ",
                   sQuote(paste("demo(package =",
                                ".packages(all.available = TRUE))")),
                   "\n",
                   "to list the demos in all *available* packages.")
	else
	    NULL
	y <- list(title = "Demos", header = NULL, results = db,
		  footer = footer)
	class(y) <- "packageIQR"
	return(y)
    }

    if(!character.only) {
    	topic <- substitute(topic)
    	if (is.call(topic) && (topic[[1L]] == "::" || topic[[1L]] == ":::")) {
	    package <- as.character(topic[[2L]])
	    topic <- as.character(topic[[3L]])
	} else
	    topic <- as.character(topic)
    }

    available <- character()
    paths <- file.path(paths, "demo")
    for(p in paths) {
	files <- basename(tools::list_files_with_type(p, "demo"))
	## Files with base names sans extension matching topic
	files <- files[topic == tools::file_path_sans_ext(files)]
	if(length(files))
	    available <- c(available, file.path(p, files))
    }
    if(length(available) == 0L)
	stop(gettextf("No demo found for topic %s", sQuote(topic)), domain = NA)
    if(length(available) > 1L) {
	available <- available[1L]
	warning(gettextf("Demo for topic %s' found more than once,\nusing the one found in %s",
                sQuote(topic), sQuote(dirname(available[1L]))), domain = NA)
    }

    ## now figure out if the package has an encoding
    pkgpath <- dirname(dirname(available))
    if (file.exists(file <- file.path(pkgpath, "Meta", "package.rds"))) {
        desc <- readRDS(file)$DESCRIPTION
        if (length(desc) == 1L) {
            enc <- as.list(desc)[["Encoding"]]
            !if(!is.null(enc)) encoding <- enc
        }
    }

    if(ask == "default")
        ask <- echo && grDevices::dev.interactive(orNone = TRUE)

    if(.Device != "null device") {
	oldask <- grDevices::devAskNewPage(ask = ask)
        on.exit(grDevices::devAskNewPage(oldask), add = TRUE)
    }

    op <- options(device.ask.default = ask)
    on.exit(options(op), add = TRUE)

    if (echo) {
	cat("\n\n",
	    "\tdemo(", topic, ")\n",
	    "\t---- ", rep.int("~", nchar(topic, type = "w")), "\n",
	    sep = "")
	if(ask && interactive())
	    readline("\nType  <Return>	 to start : ")
    }
    source(available, echo = echo, max.deparse.length = Inf,
           keep.source = TRUE, encoding = encoding)
}
#  File src/library/utils/R/edit.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

dataentry <- function (data, modes)
{
    if(!is.list(data) || !length(data) || !all(sapply(data, is.vector)))
        stop("invalid 'data' argument")
    if(!is.list(modes) ||
       (length(modes) && !all(sapply(modes, is.character))))
        stop("invalid 'modes' argument")
    .External2(C_dataentry, data, modes)
}

View <- function (x, title)
{
    ## could multi-line deparse with maliciously-designed inputs
    if(missing(title)) title <- paste("Data:", deparse(substitute(x))[1])
    as.num.or.char <- function(x)
    {
        if (is.character(x)) x
        else if (is.numeric(x)) {storage.mode(x) <- "double"; x}
        else as.character(x)
    }
    x0 <- as.data.frame(x)
    x <- lapply(x0, as.num.or.char)
    rn <- row.names(x0)
    if(any(rn != seq_along(rn))) x <- c(list(row.names = rn), x)
    if(!is.list(x) || !length(x) || !all(sapply(x, is.atomic)) ||
       !max(sapply(x, length)))
        stop("invalid 'x' argument")
    invisible(.External2(C_dataviewer, x, title))
}

edit <- function(name,...)UseMethod("edit")

edit.default <-
    function (name = NULL, file = "", title = NULL,
              editor = getOption("editor"), ...)
{
    if (is.null(title)) title <- deparse(substitute(name))
    if (is.function(editor)) invisible(editor(name, file, title))
    else .External2(C_edit, name, file, title, editor)
}

edit.data.frame <-
    function(name, factor.mode = c("character", "numeric"),
             edit.row.names =  any(row.names(name) != 1L:nrow(name)), ...)
{
    if (.Platform$OS.type == "unix"  && .Platform$GUI != "AQUA")
        if(.Platform$GUI == "unknown" || Sys.getenv("DISPLAY") == "" )
            return (edit.default(name, ...))

    is.vector.unclass <- function(x) is.vector(unclass(x))
    if (length(name) && !all(sapply(name, is.vector.unclass)
                                 | sapply(name, is.factor)))
        stop("can only handle vector and factor elements")

    factor.mode <- match.arg(factor.mode)

    as.num.or.char <- function(x)
    {
        if (is.numeric(x)) x
        else if (is.factor(x) && factor.mode == "numeric") as.numeric(x)
        else as.character(x)
    }

    attrlist <- lapply(name, attributes)
    datalist <- lapply(name, as.num.or.char)
    factors <- if (length(name))
        which(sapply(name, is.factor))
    else
        numeric()

    logicals <- if (length(name))
    	which(sapply(name, is.logical))
    else
    	numeric()

    if(length(name)) {
        has_class <-
            sapply(name, function(x) (is.object(x) || isS4(x)) && !is.factor(x))
        if(any(has_class))
            warning(sprintf(ngettext(sum(has_class),
                                    "class discarded from column %s",
                                    "classes discarded from columns %s"),
                            paste(sQuote(names(name)[has_class]),
                                  collapse=", ")),
                    domain = NA, call. = FALSE, immediate. = TRUE)
    }

    modes <- lapply(datalist, mode)
    if (edit.row.names) {
        datalist <- c(list(row.names = row.names(name)), datalist)
        modes <- c(list(row.names = "character"), modes)
    }
    rn <- attr(name, "row.names")

    out <- .External2(C_dataentry, datalist, modes)
    if(length(out) == 0L) {
        ## e.g. started with 0-col data frame or NULL, and created no cols
        return (name)
    }
    lengths <- sapply(out, length)
    maxlength <- max(lengths)
    if (edit.row.names) rn <- out[[1L]]
    for (i in which(lengths != maxlength))
         out[[i]] <- c(out[[i]], rep.int(NA, maxlength - lengths[i]))
    if (edit.row.names) {
        out <- out[-1L]
        if((ln <- length(rn)) < maxlength)
            rn <- c(rn, paste0("row", (ln+1):maxlength))
    } else if(length(rn) != maxlength) rn <- seq_len(maxlength)
    for (i in factors) {
        if(factor.mode != mode(out[[i]])) next # user might have switched mode
        a <- attrlist[[i]]
        if (factor.mode == "numeric") {
            o <- as.integer(out[[i]])
            ok <- is.na(o) | (o > 0 & o <= length(a$levels))
            if (any(!ok)) {
                warning(gettextf("invalid factor levels in '%s'", names(out)[i]),
                        domain = NA)
                o[!ok] <- NA
            }
	    attributes(o) <- a
        } else {
            o <- out[[i]]
            if (any(new <- is.na(match(o, c(a$levels, NA_integer_))))) {
                new <- unique(o[new])
                warning(gettextf("added factor levels in '%s'", names(out)[i]),
                        domain = NA)
                o <- factor(o, levels=c(a$levels, new),
                            ordered = is.ordered(o))
            } else {
                o <- match(o, a$levels)
                attributes(o) <- a
            }
        }
        out[[i]] <- o
    }
    for (i in logicals) out[[i]] <- as.logical(out[[i]])

    attr(out, "row.names") <- rn
    attr(out, "class") <- "data.frame"
    if (edit.row.names) {
        if(anyDuplicated(rn)) {
            warning("edited row names contain duplicates and will be ignored")
            attr(out, "row.names") <- seq_len(maxlength)
        }
    }
    out
}

edit.matrix <-
    function(name, edit.row.names = !is.null(dn[[1L]]), ...)
{
    if (.Platform$OS.type == "unix" && .Platform$GUI != "AQUA")
        if(.Platform$GUI == "unknown" || Sys.getenv("DISPLAY")=="" )
            return (edit.default(name, ...))
    if(!is.matrix(name) ||
       ! mode(name) %in% c("numeric", "character", "logical") ||
       any(dim(name) < 1))
        stop("invalid input matrix")
    ## logical matrices will be edited as character
    logicals <- is.logical(name)
    if (logicals) mode(name) <- "character"
    if(is.object(name) || isS4(name))
        warning("class of 'name' will be discarded",
                call. = FALSE, immediate. = TRUE)

    dn <- dimnames(name)
    datalist <- split(name, col(name))
    if(!is.null(dn[[2L]])) names(datalist) <- dn[[2L]]
    else names(datalist) <- paste0("col", 1L:ncol(name))
    modes <- as.list(rep.int(mode(name), ncol(name)))
    ## guard aginst user error (PR#10500)
    if(edit.row.names && is.null(dn[[1L]]))
        stop("cannot edit NULL row names")
    if (edit.row.names) {
        datalist <- c(list(row.names = dn[[1L]]), datalist)
        modes <- c(list(row.names = "character"), modes)
    }

    out <- .External2(C_dataentry, datalist, modes)

    lengths <- sapply(out, length)
    maxlength <- max(lengths)
    if (edit.row.names) rn <- out[[1L]]
    for (i in which(lengths != maxlength))
         out[[i]] <- c(out[[i]], rep.int(NA, maxlength - lengths[i]))
    if (edit.row.names) {
        out <- out[-1L]
        if((ln <- length(rn)) < maxlength)
            rn <- c(rn, paste0("row", (ln+1L):maxlength))
    }
    out <- do.call("cbind", out)
    if (edit.row.names)
        rownames(out) <- rn
    else if(!is.null(dn[[1L]]) && length(dn[[1L]]) == maxlength)
        rownames(out) <- dn[[1L]]
    if (logicals) mode(out) <- "logical"
    out
}

file.edit <-
  function (..., title = file, editor=getOption("editor"), fileEncoding="")
{
    file <- path.expand(c(...))
    title <- rep(as.character(title), len=length(file))
    if(nzchar(fileEncoding) && fileEncoding != "native.enc") {
        tfile <- file
        for(i in seq_along(file)) {
            ## We won't know when that is done with
            ## so leave around for the R session.
            tfile <- tempfile()
            con <- file(file[i], encoding = fileEncoding)
            writeLines(readLines(con), tfile)
            close(con)
            file[i] <- tfile
        }
    }
    if (is.function(editor)) invisible(editor(file = file, title = title))
    else invisible(.External2(C_fileedit, file, title, editor))
}

vi <- function(name = NULL, file = "")
    edit.default(name, file, editor = "vi")

emacs <- function(name = NULL, file = "")
    edit.default(name, file, editor = "emacs")

xemacs <- function(name = NULL, file = "")
    edit.default(name, file, editor = "xemacs")

xedit <- function(name = NULL, file = "")
    edit.default(name, file, editor = "xedit")

pico <- function(name = NULL, file = "")
    edit.default(name, file, editor = "pico")

#  File src/library/utils/R/example.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

## Examples as from 2.11.0 will always be new-style and hence in UTF-8
example <-
function(topic, package = NULL, lib.loc = NULL,
         character.only = FALSE, give.lines = FALSE, local = FALSE,
	 echo = TRUE, verbose = getOption("verbose"), setRNG = FALSE,
         ask = getOption("example.ask"),
	 prompt.prefix = abbreviate(topic, 6),
	 run.dontrun = FALSE)
{
    if (!character.only) {
        topic <- substitute(topic)
        if(!is.character(topic)) topic <- deparse(topic)[1L]
    }
    pkgpaths <- find.package(package, lib.loc, verbose = verbose)
    ## will only return at most one path
    file <- index.search(topic, pkgpaths, TRUE)
    if(!length(file)) {
	warning(gettextf("no help found for %s", sQuote(topic)), domain = NA)
	return(invisible())
    }
    packagePath <- dirname(dirname(file))
    pkgname <- basename(packagePath)
    lib <- dirname(packagePath)
    tf <- tempfile("Rex")
    tools::Rd2ex(.getHelpFile(file), tf, commentDontrun = !run.dontrun)
    if (!file.exists(tf)) {
	if(give.lines) return(character())
        warning(gettextf("%s has a help file but no examples", sQuote(topic)),
                domain = NA)
        return(invisible())
    }
    on.exit(unlink(tf))
    if(give.lines)
	return(readLines(tf))
    if(pkgname != "base")
        library(pkgname, lib.loc = lib, character.only = TRUE)
    if(!is.logical(setRNG) || setRNG) {
	## save current RNG state:
	if((exists(".Random.seed", envir = .GlobalEnv))) {
	    oldSeed <- get(".Random.seed", envir = .GlobalEnv)
	    on.exit(assign(".Random.seed", oldSeed, envir = .GlobalEnv),
                    add = TRUE)
	} else {
	    oldRNG <- RNGkind()
	    on.exit(RNGkind(oldRNG[1L], oldRNG[2L]), add = TRUE)
	}
	## set RNG
	if(is.logical(setRNG)) { # i.e. == TRUE: use the same as R CMD check
	    ## see share/R/examples-header.R
	    RNGkind("default", "default")
	    set.seed(1)
	} else eval(setRNG)
    }
    zz <- readLines(tf, n = 1L)
    skips <- 0L
    if (echo) {
	## skip over header
	zcon <- file(tf, open="rt")
	while(length(zz) && !length(grep("^### \\*\\*", zz))) {
	    skips <- skips + 1L
	    zz <- readLines(zcon, n=1L)
	}
	close(zcon)
    }
    if(ask == "default")
        ask <- echo && grDevices::dev.interactive(orNone = TRUE)
    if(ask) {
	if(.Device != "null device") {
	    oldask <- grDevices::devAskNewPage(ask = TRUE)
            if(!oldask) on.exit(grDevices::devAskNewPage(oldask), add = TRUE)
        }
        ## <FIXME>
        ## This ensures that any device opened by the examples will
        ## have ask = TRUE set, but it does not return the device to
        ## the expected 'ask' state if it is left as the current device.
        ## </FIXME>
        op <- options(device.ask.default = TRUE)
        on.exit(options(op), add = TRUE)
    }
    source(tf, local, echo = echo,
           prompt.echo = paste0(prompt.prefix, getOption("prompt")),
           continue.echo = paste0(prompt.prefix, getOption("continue")),
           verbose = verbose, max.deparse.length = Inf, encoding = "UTF-8",
    	   skip.echo = skips, keep.source=TRUE)
}
#  File src/library/utils/R/filetest.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

### ** file_test

file_test <-
function(op, x, y)
{
    ## Provide shell-style '-f', '-d', '-x', '-nt' and '-ot' tests.
    ## Note that file.exists() only tests existence ('test -e' on some
    ## systems), and that our '-f' tests for existence and not being a
    ## directory (the GNU variant tests for being a regular file).
    ## Note: vectorized in x and y.
    switch(op,
           "-f" = !is.na(isdir <- file.info(x)$isdir) & !isdir,
           "-d" = !is.na(isdir <- file.info(x)$isdir) & isdir,
           "-nt" = (!is.na(mt.x <- file.info(x)$mtime)
                    & !is.na(mt.y <- file.info(y)$mtime)
                    & (mt.x > mt.y)),
           "-ot" = (!is.na(mt.x <- file.info(x)$mtime)
                    & !is.na(mt.y <- file.info(y)$mtime)
                    & (mt.x < mt.y)),
           "-x" = (file.access(x, 1L) == 0L),
           stop(gettextf("test '%s' is not available", op),
                domain = NA))
}
#  File src/library/utils/R/fineLineNum.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 2009-2012 Duncan Murdoch and the R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


.normalizePath <- function(path, wd) {
    if (!missing(wd) && !is.null(wd)) {
    	oldwd <- setwd(wd)
    	on.exit(setwd(oldwd))
    }
    suppressWarnings(normalizePath(path))
}

fnLineNum <- function(f, srcfile, line, nameonly=TRUE) {

    stopifnot(length(line) == 1)

    targetfilename <- .normalizePath(srcfile$filename)

    fnsrc <- attr(f, "srcref")
    if (!is.null(fnsrc))
    	fnsrc <- attr(fnsrc, "srcfile")
    else
    	fnsrc <- attr(body(f), "srcfile")
    if (is.null(fnsrc)) return(NULL)

    if (missing(srcfile)) {
    	srcfile <- fnsrc
    }

    isBrace <- function(expr)
        typeof(expr) == "symbol" && identical(as.character(expr), "{")

    lineNumInExpr <- function(expr, haveSrcrefs = FALSE) {
	if (typeof(expr) == "language") {
	    srcrefs <- attr(expr, "srcref")
	    for (i in seq_along(expr)) {
		srcref <- srcrefs[[i]]
		# Check for non-matching range
		if (!is.null(srcref) && (srcref[1] > line || line > srcref[3]))  next
		# We're in range.  See if there's a finer division
		finer <- lineNumInExpr(expr[[i]], haveSrcrefs || !is.null(srcrefs))
		if (!is.null(finer)) {
		    return(c(i, finer))
		}
		# Do we have a srcref?  It must point to this expression.
		# But do avoid matching the opening brace in a block:  match the whole block
		# instead.

		havebrace <- isBrace(expr[[i]])
		if (!is.null(srcref)
		    && (!haveSrcrefs || !havebrace)) {
		    return(i)
		}
	    }
	}
	return(NULL)
    }

    perfectMatch <- identical(.normalizePath(fnsrc$filename, fnsrc$wd), targetfilename)
    if (perfectMatch ||
        (nameonly && !is.null(fnsrc$filename) && basename(fnsrc$filename) == basename(targetfilename))) {
	if (!is.na(srcfile$timestamp) && fnsrc$timestamp != srcfile$timestamp)
	    timediff <- fnsrc$timestamp - srcfile$timestamp
	else
	    timediff <- 0
	at <- lineNumInExpr(body(f))
	if (!is.null(at))
	  return(list(at=at, filename=.normalizePath(fnsrc$filename, fnsrc$wd), line=line,
	              timediff=timediff))
    }
    return(NULL)
}

findLineNum <- function(srcfile, line, nameonly=TRUE, envir=parent.frame(),
			lastenv) {
    count <- 0
    result <- list()

    if (!inherits(srcfile, "srcfile")) {
    	if (missing(line)) {
    	    line <- as.numeric(sub(".*#", "", srcfile))
    	    if (is.na(line)) stop("Line number missing")
    	    srcfile <- sub("#[^#]*", "", srcfile)
    	}

    	srcfile <- srcfile(srcfile)
    }

    if (missing(lastenv)) {
    	if (missing(envir)) lastenv <- globalenv()
    	else lastenv <- emptyenv()
    }

    if (!is.environment(envir))
    	envir <- environment(envir)

    fns <- character()
    envirs <- list()
    e <- envir
    repeat {
    	fns <- c(fns, lsf.str(envir=e, all=TRUE))
    	oldlen <- length(envirs)
    	length(envirs) <- length(fns)
    	if (length(envirs) > oldlen)
    	    for (i in seq.int(oldlen+1, length(envirs))) envirs[[i]] <- e
    	if (identical(e, lastenv) || identical(e, emptyenv())) break
    	e <- parent.env(e)
    }

    for (i in seq_along(fns)) {
	functionName <- fns[i]
	fn <- get(functionName, envir=envirs[[i]])
	loc <- fnLineNum(fn, srcfile=srcfile, line=line,
    	                  nameonly=nameonly)
    	if (!is.null(loc)) {
    	    count <- count + 1
    	    result[[count]] <- c(list(name=functionName, env=envirs[[i]]), loc)
    	}
    	gen <- tryCatch(methods::isGeneric(functionName, envirs[[i]], fdef=fn),
                        error = identity)
    	if (isTRUE(gen)) {
    	    e1 <- environment(fn)$.AllMTable
    	    if (!is.null(e1)) {
		sigs <- ls(e1)
		for (j in seq_along(sigs)) {
		    sig <- sigs[j]
		    fn <- get(sig, e1)
		    if (typeof(fn) != "closure") next

		    loc <- fnLineNum(fn, srcfile=srcfile, line=line,
				    nameonly=nameonly)
		    if (is.null(loc)
		        && length(body(fn)) > 1
		        && length(body(fn)[[2]]) > 2
		        && typeof(body(fn)[[c(2,3)]]) == "closure") {
				# desperate try:  look for
		    		# .local <- original defn
		    	fn2 <- body(fn)[[c(2,3)]]

		    	loc <- fnLineNum(fn2, srcfile=srcfile, line=line,
				    nameonly=nameonly)
		 	# FIXME:  can trace() set a breakpoint
		 	#	  within a function like this?
		        if (!is.null(loc)) loc$at <- c(2,3)
		    }
		    if (!is.null(loc)) {
			count <- count + 1
			result[[count]] <- c(list(name=functionName, env=envirs[[i]],
						signature=strsplit(sig, "#")[[1]]), loc)
		    }
		}
	    }
	}
    }
    return(structure(result, class="findLineNumResult"))
}

print.findLineNumResult <- function(x, steps=TRUE, ...) {
    if (!length(x)) cat("No source refs found.\n")
    filename <- NULL
    line <- 0
    for (i in seq_along(x)) {
    	if (!identical(filename, x[[i]]$filename) ||
    	    !identical(line, x[[i]]$line)) {
    	    filename <- x[[i]]$filename
    	    line <- x[[i]]$line
    	    cat(filename, "#", line, ":\n", sep = "")
    	}
        cat(" ", x[[i]]$name, if (steps) paste(" step ", paste(x[[i]]$at, collapse=",")) else "", sep = "")
        if (!is.null(x[[i]]$signature))
            cat(" signature ", paste(x[[i]]$signature, collapse=","), sep = "")
        cat(" in ", format(x[[i]]$env), "\n", sep = "")
    }
}


setBreakpoint <- function(srcfile, line, nameonly=TRUE, envir=parent.frame(), lastenv,
                          verbose = TRUE, tracer, print=FALSE, clear=FALSE,
                         ...) {

    if (missing(lastenv)) {
    	if (missing(envir)) lastenv <- globalenv()
    	else lastenv <- emptyenv()
    }
    locations <- findLineNum(srcfile, line, nameonly, envir, lastenv)
    if (verbose) print(locations, steps=!clear)
    breakpoint <- missing(tracer)
    while (length(locations)) {
    	what <- locations[[1]]$name
    	where <- locations[[1]]$env
    	at <- list(locations[[1]]$at)
    	signature <- locations[[1]]$signature
    	if (breakpoint) {
    	    filename <- basename(locations[[1]]$filename)
    	    linenum <- locations[[1]]$line
	    tracer <- bquote({cat(paste0(.(filename), "#", .(linenum), "\n"))
    	                      browser(skipCalls=4L)})
    	}
    	locations[[1]] <- NULL
        i <- 1
    	while (i <= length(locations)) {
    	    if (what == locations[[i]]$name &&
    	        identical(where, locations[[i]]$env) &&
    	        identical(signature, locations[[i]]$signature))	 {
    	    	at <- c(at, list(locations[[i]]))
    	    	locations[[i]] <- NULL
    	    } else
    	    	i <- i+1
    	}
    	if (clear) {
    	    if (is.null(signature))
  		untrace(what, where=where)
    	    else
    	    	untrace(what, signature=signature, where=where)
    	} else if (is.null(signature))
    	    trace(what, tracer, at=at, where=where, print=print, ...)
    	else
    	    trace(what, signature=signature, tracer, at=at, where=where, ...)
    }
}
#  File src/library/utils/R/fix.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

fix <- function (x, ...)
{
    subx <- substitute(x)
    if (is.name(subx))
        subx <- deparse(subx)
    if (!is.character(subx) || length(subx) != 1L)
        stop("'fix' requires a name")
    parent <- parent.frame()
    if (exists(subx, envir=parent, inherits = TRUE))
        x <- edit(get(subx, envir=parent), title = subx, ...)
    else {
        x <- edit(function(){}, title = subx, ...)
        environment(x) <- .GlobalEnv
    }
    assign(subx, x, envir = .GlobalEnv)
}
#  File src/library/utils/R/format.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

formatUL <-
function(x, label = "*", offset = 0,
         width = 0.9 * getOption("width"))
{
    if(!length(x))
        return(character())
    .format_rl_table(label, x, offset, width)
}

formatOL <-
function(x, type = "arabic", offset = 0, start = 1,
         width = 0.9 * getOption("width"))
{
    if(!length(x))
        return(character())
    type_tokens <- c("1", "A", "a", "I", "i")
    type_full_names <- c("arabic", "Alph", "alph", "Roman", "roman")
    type <- match.arg(type, c(type_tokens, type_full_names))
    if(nchar(type, "b") > 1L)
        type <- type_tokens[match(type, type_full_names)]
    len <- length(x)
    labels <- seq.int(start[1L], length.out = len)
    upper <- labels[len]
    if(type %in% c("A", "a")) {
        if(upper > 26L)
            stop(gettextf("too many list items (at most up to %d)", 26L),
                 domain = NA)
        labels <- if(type == "A")
            LETTERS[labels]
        else
            letters[labels]
    }
    else if(type %in% c("I", "i")) {
        if(upper > 3899L)
            stop(gettextf("too many list items (at most up to %d)", 3899L),
                 domain = NA)
        labels <- as.character(as.roman(labels))
        if(type == "i")
            labels <- tolower(labels)
    }
    .format_rl_table(sprintf("%s.", labels), x, offset, width)
}

.format_rl_table <-
function(labels, x, offset = 0, width = 0.9 * getOption("width"),
         sep = " ")
{
    ## Format a 2-column table with right-justified item labels and
    ## left-justified text.  Somewhat tricky because strwrap() eats up
    ## leading whitespace ...

    .make_empty_string <- function(n) {
        paste(rep.int(" ", n), collapse = "")
    }

    labels <- format(labels, justify = "right")
    len <- length(x)
    delta <- nchar(labels[1L], "width") + offset
    x <- strwrap(x, width = width - delta - nchar(sep, "width"),
                 simplify = FALSE)
    nlines <- cumsum(sapply(x, length))
    prefix <- rep.int(.make_empty_string(delta), nlines[len])
    prefix[1L + c(0L, nlines[-len])] <-
        paste0(.make_empty_string(offset), labels)
    paste(prefix, unlist(x), sep = sep)
}
#  File src/library/utils/R/frametools.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

stack <- function(x, ...) UseMethod("stack")

stack.data.frame <- function(x, select, ...)
{
    if (!missing(select)) {
	nl <- as.list(1L:ncol(x))
	names(nl) <- names(x)
	vars <- eval(substitute(select),nl, parent.frame())
        x <- x[, vars, drop=FALSE]
    }
    keep <- unlist(lapply(x, is.vector))
    if(!sum(keep)) stop("no vector columns were selected")
    if(!all(keep))
        warning("non-vector columns will be ignored")
    x <- x[, keep, drop = FALSE]
    ## need to avoid promotion to factors
    data.frame(values = unlist(unname(x)),
               ind = factor(rep.int(names(x), lapply(x, length))),
               stringsAsFactors = FALSE)
}

stack.default <- function(x, ...)
{
    x <- as.list(x)
    keep <- unlist(lapply(x, is.vector))
    if(!sum(keep)) stop("at least one vector element is required")
    if(!all(keep)) warning("non-vector elements will be ignored")
    x <- x[keep]
    data.frame(values = unlist(unname(x)),
               ind = factor(rep.int(names(x), lapply(x, length))),
               stringsAsFactors = FALSE)
}

unstack <- function(x, ...) UseMethod("unstack")

unstack.data.frame <- function(x, form, ...)
{
    form <- if(missing(form)) stats::formula(x) else stats::as.formula(form)
    if (length(form) < 3)
        stop("'form' must be a two-sided formula")
    res <- c(tapply(eval(form[[2L]], x), eval(form[[3L]], x), as.vector))
    if (length(res) >= 2L && any(diff(unlist(lapply(res, length))) != 0L))
        return(res)
    data.frame(res, stringsAsFactors = FALSE)
}

unstack.default <- function(x, form, ...)
{
    x <- as.list(x)
    form <- stats::as.formula(form)
    if ((length(form) < 3) || (length(all.vars(form))>2))
        stop("'form' must be a two-sided formula with one term on each side")
    res <- c(tapply(eval(form[[2L]], x), eval(form[[3L]], x), as.vector))
    if (length(res) >= 2L && any(diff(unlist(lapply(res, length))) != 0L))
        return(res)
    data.frame(res, stringsAsFactors = FALSE)
}
#  File src/library/utils/R/glob2rx.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

glob2rx <- function(pattern, trim.head = FALSE, trim.tail = TRUE)
{
    ## Purpose: Change 'ls' aka 'wildcard' aka 'globbing' _pattern_ to
    ##	      Regular Expression (as in grep, perl, emacs, ...)
    ## -------------------------------------------------------------------------
    ## Author: Martin Maechler ETH Zurich, ~ 1991
    ##	       New version using [g]sub() : 2004
    p <- gsub("\\.","\\\\.", paste0("^", pattern, "$"))
    p <- gsub("\\?",	 ".",  gsub("\\*",  ".*", p))
    ## 'Escaping hell' : at least for '(', '[' and '{'
    p <- gsub("([^\\])\\(", "\\1\\\\(", p)
    p <- gsub("([^\\])\\[", "\\1\\\\[", p)
    p <- gsub("([^\\])\\{", "\\1\\\\{", p)
    ## these are trimming ".*$" and "^.*" - in most cases only for aesthetics
    if(trim.tail) p <- sub("\\.\\*\\$$", "", p)
    if(trim.head) p <- sub("\\^\\.\\*",  "", p)
    p
}
#  File src/library/utils/R/head.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

### placed in the public domain 2002
### Patrick Burns patrick@burns-stat.com
###
### Adapted for negative arguments by Vincent Goulet
### <vincent.goulet@act.ulaval.ca>, 2006

head <- function(x, ...) UseMethod("head")

head.default <- function(x, n = 6L, ...)
{
    stopifnot(length(n) == 1L)
    n <- if (n < 0L) max(length(x) + n, 0L) else min(n, length(x))
    x[seq_len(n)]
}

## head.matrix and tail.matrix are now exported (to be used for other classes)
head.data.frame <- head.matrix <- function(x, n = 6L, ...)
{
    stopifnot(length(n) == 1L)
    n <- if (n < 0L) max(nrow(x) + n, 0L) else min(n, nrow(x))
    x[seq_len(n), , drop=FALSE]
}
head.table  <- function(x, n = 6L, ...) {
    (if(length(dim(x)) == 2L) head.matrix else head.default)(x, n=n)
}

head.ftable <- function(x, n = 6L, ...) {
    r <- format(x)
    dimnames(r) <- list(rep.int("", nrow(r)), rep.int("", ncol(r)))
    noquote(head.matrix(r, n = n + nrow(r) - nrow(x), ...))
}

head.function <- function(x, n = 6L, ...)
{
    lines <- as.matrix(deparse(x))
    dimnames(lines) <- list(seq_along(lines),"")
    noquote(head(lines, n=n))
}

tail <- function(x, ...) UseMethod("tail")

tail.default <- function(x, n = 6L, ...)
{
    stopifnot(length(n) == 1L)
    xlen <- length(x)
    n <- if (n < 0L) max(xlen + n, 0L) else min(n, xlen)
    x[seq.int(to = xlen, length.out = n)]
}

tail.data.frame <- function(x, n = 6L, ...)
{
    stopifnot(length(n) == 1L)
    nrx <- nrow(x)
    n <- if (n < 0L) max(nrx + n, 0L) else min(n, nrx)
    x[seq.int(to = nrx, length.out = n), , drop = FALSE]
}

tail.matrix <- function(x, n = 6L, addrownums = TRUE, ...)
{
    stopifnot(length(n) == 1L)
    nrx <- nrow(x)
    n <- if (n < 0L) max(nrx + n, 0L) else min(n, nrx)
    sel <- seq.int(to = nrx, length.out = n)
    ans <- x[sel, , drop = FALSE]
    if (addrownums && is.null(rownames(x)))
    	rownames(ans) <- paste0("[", sel, ",]")
    ans
}
tail.table  <- function(x, n = 6L, addrownums = TRUE, ...) {
    (if(length(dim(x)) == 2L) tail.matrix else tail.default)(x, n=n,
	      addrownums = addrownums, ...)
}

tail.ftable <- function(x, n = 6L, addrownums = FALSE, ...) {
    r <- format(x)
    dimnames(r) <- list(if(!addrownums) rep.int("", nrow(r)),
			rep.int("", ncol(r)))
    noquote(tail.matrix(r, n = n, addrownums = addrownums, ...))
}

tail.function <- function(x, n = 6L, ...)
{
    lines <- as.matrix(deparse(x))
    dimnames(lines) <- list(seq_along(lines),"")
    noquote(tail(lines, n=n))
}
#  File src/library/utils/R/help.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

help <-
function(topic, package = NULL, lib.loc = NULL,
         verbose = getOption("verbose"),
         try.all.packages = getOption("help.try.all.packages"),
         help_type = getOption("help_type"))
{
    types <- c("text", "html", "pdf")
    if(!missing(package))
        if(is.name(y <- substitute(package)))
            package <- as.character(y)

    ## If no topic was given ...
    if(missing(topic)) {
        if(!missing(package)) {         # "Help" on package.
            help_type <- if(!length(help_type)) "text"
            else match.arg(tolower(help_type), types)
            ## Carter Butts and others misuse 'help(package=)' in startup
            if (interactive() && help_type == "html") {
                if (tools:::httpdPort == 0L) tools::startDynamicHelp()
                if (tools:::httpdPort <= 0L) # fallback to text help
                    return(library(help = package, lib.loc = lib.loc,
                                   character.only = TRUE))
                browser <- if (.Platform$GUI == "AQUA") {
                    get("aqua.browser", envir = as.environment("tools:RGUI"))
                } else getOption("browser")
 		browseURL(paste0("http://127.0.0.1:", tools:::httpdPort,
                                 "/library/", package, "/html/00Index.html"),
                          browser)
                return(invisible())
            } else return(library(help = package, lib.loc = lib.loc,
                                  character.only = TRUE))
        }
        if(!missing(lib.loc))           # text "Help" on library.
            return(library(lib.loc = lib.loc))
        ## ultimate default is to give help on help()
        topic <- "help"; package <- "utils"; lib.loc <- .Library
    }

    ischar <- tryCatch(is.character(topic) && length(topic) == 1L,
                       error = identity)
    if(inherits(ischar, "error")) ischar <- FALSE
    ## if this was not a length-one character vector, try for the name.
    if(!ischar) {
        ## the reserved words that could be parsed as a help arg:
        reserved <-
            c("TRUE", "FALSE", "NULL", "Inf", "NaN", "NA", "NA_integer_",
              "NA_real_", "NA_complex_", "NA_character_")
        stopic <- deparse(substitute(topic))
        if(!is.name(substitute(topic)) && ! stopic %in% reserved)
            stop("'topic' should be a name, length-one character vector or reserved word")
        topic <- stopic
    }

    help_type <- if(!length(help_type)) "text"
    else match.arg(tolower(help_type), types)

    paths <- index.search(topic,
                          find.package(package, lib.loc, verbose = verbose))
    tried_all_packages <- FALSE
    if(!length(paths)
       && is.logical(try.all.packages) && !is.na(try.all.packages)
       && try.all.packages && missing(package) && missing(lib.loc)) {
        ## Try all the remaining packages.
        for(lib in .libPaths()) {
            packages <- .packages(TRUE, lib)
            packages <- packages[is.na(match(packages, .packages()))]
            paths <- c(paths, index.search(topic, file.path(lib, packages)))
        }
        paths <- paths[paths != ""]
        tried_all_packages <- TRUE
    }

    paths <- unique(paths)
    attributes(paths) <-
        list(call = match.call(), topic = topic,
             tried_all_packages = tried_all_packages, type = help_type)
    class(paths) <- "help_files_with_topic"
    paths
}

print.help_files_with_topic <- function(x, ...)
{
    browser <- getOption("browser")
    topic <- attr(x, "topic")
    type <- attr(x, "type")
    if (.Platform$GUI == "AQUA" && type == "html")
        browser <- get("aqua.browser", envir = as.environment("tools:RGUI"))
    paths <- as.character(x)
    if(!length(paths)) {
        writeLines(c(gettextf("No documentation for %s in specified packages and libraries:",
                              sQuote(topic)),
                     gettextf("you could try %s",
                              sQuote(paste0("??", topic)))))
        return(invisible(x))
    }

    if(type == "html")
        if (tools:::httpdPort == 0L) tools::startDynamicHelp()

    if(attr(x, "tried_all_packages")) {
        paths <- unique(dirname(dirname(paths)))
        msg <- gettextf("Help for topic %s is not in any loaded package but can be found in the following packages:",
                        sQuote(topic))
        if (type == "html" && tools:::httpdPort > 0L) {
            path <- file.path(tempdir(), ".R/doc/html")
            dir.create(path, recursive = TRUE, showWarnings = FALSE)
            out <- paste0('<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n',
                          '<html><head><title>R: help</title>\n',
                          '<meta http-equiv="Content-Type" content="text/html; charset="UTF-8">\n',
                          '<link rel="stylesheet" type="text/css" href="/doc/html/R.css">\n',
                          '</head><body>\n\n<hr>\n')
            out <- c(out, '<p>', msg, '</p><br>')
            out <- c(out, '<table width="100%" summary="R Package list">\n',
                     '<tr align="left" valign="top">\n',
                     '<td width="25%">Package</td><td>Library</td></tr>\n')
            pkgs <- basename(paths)
            links <- paste0('<a href="http://127.0.0.1:', tools:::httpdPort,
                            '/library/', pkgs, '/help/', topic, '">',
                            pkgs, '</a>')
            out <- c(out, paste0('<tr align="left" valign="top">\n',
                                '<td>', links, '</td><td>',
                                dirname(paths), '</td></tr>\n'))
            out <- c(out, "</table>\n</p>\n<hr>\n</body></html>")
            writeLines(out, file.path(path, "all.available.html"))
            browseURL(paste0("http://127.0.0.1:", tools:::httpdPort,
                             "/doc/html/all.available.html"), browser)
        } else {
            writeLines(c(strwrap(msg), "",
                         paste(" ",
                               formatDL(c(gettext("Package"), basename(paths)),
                                        c(gettext("Library"), dirname(paths)),
                                        indent = 22))))
        }
    } else {
        if(length(paths) > 1L) {
            if (type == "html" && tools:::httpdPort > 0L) { # Redo the search if dynamic help is running
		browseURL(paste0("http://127.0.0.1:", tools:::httpdPort,
                                 "/library/NULL/help/", topic), browser)
		return(invisible(x))
	    }
            file <- paths[1L]
            p <- paths
            msg <- gettextf("Help on topic %s was found in the following packages:",
                            sQuote(topic))
            paths <- dirname(dirname(paths))
            txt <- formatDL(c("Package", basename(paths)),
                            c("Library", dirname(paths)),
                            indent = 22L)
            writeLines(c(strwrap(msg), "", paste(" ", txt), ""))
            if(interactive()) {
                fp <- file.path(paths, "Meta", "Rd.rds")
                tp <- basename(p)
                titles <- tp
                if(type == "html" || type == "latex")
                    tp <- tools::file_path_sans_ext(tp)
                for (i in seq_along(fp)) {
                    tmp <- try(readRDS(fp[i]))
                    titles[i] <- if(inherits(tmp, "try-error"))
                        "unknown title" else
                    tmp[tools::file_path_sans_ext(tmp$File) == tp[i], "Title"]
                }
                txt <- paste0(titles, " {", basename(paths), "}")
                ## the default on menu() is currtently graphics = FALSE
                res <- menu(txt, title = gettext("Choose one"),
                            graphics = getOption("menu.graphics"))
                if(res > 0) file <- p[res]
            } else {
                writeLines(gettext("\nUsing the first match ..."))
            }
        }
        else
            file <- paths

        if(type == "html") {
            if (tools:::httpdPort > 0L) {
		path <- dirname(file)
		dirpath <- dirname(path)
		pkgname <- basename(dirpath)
		browseURL(paste0("http://127.0.0.1:", tools:::httpdPort,
                                 "/library/", pkgname, "/html/", basename(file),
                                 ".html"), browser)
            } else {
                warning("HTML help is unavailable", call. = FALSE)
                att <- attributes(x)
                xx <- sub("/html/([^/]*)\\.html$", "/help/\\1", x)
                attributes(xx) <- att
                attr(xx, "type") <- "text"
                print(xx)
            }
        } else if(type == "text") {
            pkgname <- basename(dirname(dirname(file)))
            temp <- tools::Rd2txt(.getHelpFile(file), out = tempfile("Rtxt"),
                                  package = pkgname)
            file.show(temp, title = gettextf("R Help on %s", sQuote(topic)),
                      delete.file = TRUE)
        }
        else if(type %in% "pdf") {
            path <- dirname(file)
            dirpath <- dirname(path)
            texinputs <- file.path(dirpath, "help", "figures")
            tf2 <- tempfile("Rlatex")
            tools::Rd2latex(.getHelpFile(file), out = tf2)
            .show_help_on_topic_offline(tf2, topic, type, texinputs)
            unlink(tf2)
        }
    }

    invisible(x)
}

.show_help_on_topic_offline <-
    function(file, topic, type = "pdf", texinputs = NULL)
{
    encoding <-""
    lines <- readLines(file)
    encpatt <- "^\\\\inputencoding\\{(.*)\\}$"
    if(length(res <- grep(encpatt, lines, perl = TRUE, useBytes = TRUE)))
        encoding <- sub(encpatt, "\\1", lines[res],
                        perl = TRUE, useBytes = TRUE)
    texfile <- paste0(topic, ".tex")
    on.exit(unlink(texfile)) ## ? leave to helper
    if(nzchar(opt <- Sys.getenv("R_RD4PDF"))) opt else "times,inconsolata"
    has_figure <- any(grepl("\\Figure", lines))
    cat("\\documentclass[", getOption("papersize"), "paper]{article}\n",
        "\\usepackage[", opt, "]{Rd}\n",
        if(nzchar(encoding)) sprintf("\\usepackage[%s]{inputenc}\n", encoding),
        "\\InputIfFileExists{Rhelp.cfg}{}{}\n",
        "\\usepackage{graphicx}\n",
        "\\begin{document}\n",
        file = texfile, sep = "")
    file.append(texfile, file)
    cat("\\end{document}\n", file = texfile, append = TRUE)
    helper <- if (exists("offline_help_helper", envir = .GlobalEnv))
        get("offline_help_helper", envir = .GlobalEnv)
    else utils:::offline_help_helper
    if (has_figure) helper(texfile, type, texinputs)
    else helper(texfile, type)
    invisible()
}


.getHelpFile <- function(file)
{
    path <- dirname(file)
    dirpath <- dirname(path)
    if(!file.exists(dirpath))
        stop(gettextf("invalid %s argument", sQuote("file")), domain = NA)
    pkgname <- basename(dirpath)
    RdDB <- file.path(path, pkgname)
    if(!file.exists(paste(RdDB, "rdx", sep = ".")))
        stop(gettextf("package %s exists but was not installed under R >= 2.10.0 so help cannot be accessed", sQuote(pkgname)), domain = NA)
    tools:::fetchRdDB(RdDB, basename(file))
}


offline_help_helper <- function(texfile, type, texinputs = NULL)
{
    ## Some systems have problems with texfile names like ".C.tex"
    tf <- tempfile("tex", tmpdir = ".", fileext = ".tex"); on.exit(unlink(tf))
    file.copy(texfile, tf)
    tools::texi2pdf(tf, clean = TRUE, texinputs = texinputs)
    ofile <- sub("tex$", "pdf", tf)
    ofile2 <- sub("tex$", "pdf", texfile)
    if(!file.exists(ofile))
        stop(gettextf("creation of %s failed", sQuote(ofile2)), domain = NA)
    if(file.copy(ofile, ofile2, overwrite = TRUE)) {
        unlink(ofile)
        message(gettextf("Saving help page to %s", sQuote(basename(ofile2))),
                domain = NA)
    } else {
        message(gettextf("Saving help page to %s", sQuote(ofile)), domain = NA)
    }
    invisible()
}

#  File src/library/utils/R/unix/help.request.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

help.request <- function (subject = "", address = "r-help@R-project.org",
			  file = "R.help.request", ...)
{
    no <- function(answer) answer == "n"
    yes <- function(answer) answer == "y"
    webpage <- "corresponding web page"
    catPlease <- function()
	cat("Please do this first - the",
	    webpage,"has been loaded in your web browser\n")
    go <- function(url) {
	catPlease()
	browseURL(url)
    }
    readMyLine <- function(..., .A. = "(y/n)")
	readline(paste(paste(strwrap(paste(...)), collapse="\n"),
		       .A., "")) # space after question

    checkPkgs <- function(pkgDescs,
			  pkgtxt = paste("packages",
                          paste(names(pkgDescs), collapse=", ")))
    {
        cat("Checking if", pkgtxt, "are up-to-date; may take some time...\n")

        stopifnot(sapply(pkgDescs, inherits, what="packageDescription"))
        fields <- .instPkgFields(NULL)
	n <- length(pkgDescs)
	iPkgs <- matrix(NA_character_, n, 2L + length(fields),
		      dimnames=list(NULL, c("Package", "LibPath", fields)))
	for(i in seq_len(n)) {
	    desc <- c(unlist(pkgDescs[[i]]),
		      "LibPath" = dirname(dirname(dirname(attr(pkgDescs[[i]],
		      "file")))))
	    nms <- intersect(names(desc), colnames(iPkgs))
	    iPkgs[i, nms] <- desc[nms]
	}

	old <- old.packages(instPkgs = iPkgs)

	if (!is.null(old)) {
	    update <- readMyLine("The following installed packages are out-of-date:\n",
				 paste(strwrap(rownames(old),
					       width = 0.7 *getOption("width"),
					       indent = 0.15*getOption("width")),
				       collapse="\n"),
				 "would you like to update now?")
	    if (yes(update)) update.packages(oldPkgs = old, ask = FALSE)
	}
    }

    cat("Checklist:\n")
    post <- readline("Have you read the posting guide? (y/n) ")
    if (no(post)) return(go("http://www.r-project.org/posting-guide.html"))
    FAQ <- readline("Have you checked the FAQ? (y/n) ")
    if (no(FAQ)) return(go("http://cran.r-project.org/faqs.html"))
    intro <- readline("Have you checked An Introduction to R? (y/n) ")
    if (no(intro))
	return(go("http://cran.r-project.org/manuals.html"))
    NEWS <- readMyLine("Have you checked the NEWS of the latest development release?")
    if (no(NEWS)) return(go("http://cran.r-project.org/doc/manuals/r-devel/NEWS.html"))
    rsitesearch <- readline("Have you looked on RSiteSearch? (y/n) ")
    if (no(rsitesearch)) {
	catPlease()
	return(RSiteSearch(subject))
    }
    inf <- sessionInfo()
    if ("otherPkgs" %in% names(inf)) {
	oPkgs <- names(inf$otherPkgs)
        ## FIXME: inf$otherPkgs is a list of packageDescription()s
	other <-
	    readMyLine("You have packages",
                       paste0("(", paste(sQuote(oPkgs), collapse=", "),")"),
                       "other than the base packages loaded. ",
		       "If your query relates to one of these, have you ",
		       "checked any corresponding books/manuals and",
		       "considered contacting the package maintainer?",
                       .A. = "(y/n/NA)")
	if(no(other)) return("Please do this first.")
    }

    page <- url("http://cran.r-project.org/bin/windows/base")
    title <- grep("<title>", readLines(page, 10L), fixed = TRUE, value = TRUE)
    ver <- sub("^.*R-([^ ]*) for Windows.*$", "\\1", title)
    if (getRversion() < numeric_version(ver)) {
	update <- readMyLine("Your R version is out-of-date,",
			     "would you like to update now?")
	if(yes(update)) return(go(getOption("repos")))
    }
    if ("otherPkgs" %in% names(inf)) {
        checkPkgs(inf$otherPkgs)
    }
    ## To get long prompt!
    cat("Have you written example code that is\n",
	"- minimal\n - reproducible\n - self-contained\n - commented",
	"\nusing data that is either\n",
	"- constructed by the code\n - loaded by data()\n",
	"- reproduced using dump(\"mydata\", file = \"\")\n")
    code <- readMyLine("have you checked this code in a fresh R session",
		       "(invoking R with the --vanilla option if possible)",
		       "and is this code copied to the clipboard?")
    if (no(code))
	return(cat("\nIf your query is not directly related to code",
		   "(e.g. a general query \nabout R's capabilities),",
		   "email R-help@r-project.org directly. ",
		   "\nOtherwise prepare some example code first.\n"))
    change <- readline(paste("Would you like to change your subject line:",
			     subject, "to something more meaningful? (y/n) ",
			     sep = "\n"))
    if (yes(change))
	subject <- readline("Enter subject: \n")

    create.post(instructions = paste(
		"\\n<<SEND AS PLAIN TEXT!>>\\n\\n",
		"\\n<<Write your query here, using your example code to illustrate>>",
		"\\n<<End with your name and affiliation>>\\n\\n\\n\\n"),
		description = "help request",
		subject = subject, address = address,
                filename = file, info = bug.report.info(), ...)
}
#  File src/library/utils/R/help.search.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

.hsearch_db <- local({
    hdb <- NULL
    function(new) {
	if(!missing(new))
	    hdb <<- new
	else
	    hdb
    }
})

merge.vignette.index <- function(hDB, path, pkg) {
    ## Vignettes in the hsearch index started in R2.14.0
    ## Most packages don't have them, so the following should not be
    ## too inefficient
    if(file.exists(v_file <- file.path(path, "Meta", "vignette.rds"))
       && !is.null(vDB <- readRDS(v_file))
       && nrow(vDB)) {
	## Make it look like an hDB base matrix and append it
	base <- matrix("", nrow=nrow(vDB), ncol=8)
	colnames(base) <- colnames(hDB[[1]])
	base[,"Package"] <- pkg
	base[,"LibPath"] <- path
	id <- as.character(1:nrow(vDB) + NROW(hDB[[1]]))
	base[,"ID"] <- id
	base[,"name"] <- sub("\\.[^.]*$", "", basename(vDB$File))
	base[,"topic"] <- base[,"name"]
	base[,"title"] <- vDB$Title
	base[,"Type"] <- "vignette"
	hDB[[1L]] <- rbind(hDB[[1L]], base)
	aliases <- matrix("", nrow=nrow(vDB), ncol=3)
	colnames(aliases) <- colnames(hDB[[2]])
	aliases[,"Aliases"] <- base[,"name"]
	aliases[,"ID"] <- id
	aliases[,"Package"] <- pkg
	hDB[[2L]] <- rbind(hDB[[2L]], aliases)
	nkeywords <- sum(sapply(vDB$Keywords, length))
	if (nkeywords) {
	    keywords <- matrix("", nrow=nkeywords, ncol=3)
	    colnames(keywords) <- colnames(hDB[[4]])
	    keywords[,"Concepts"] <- unlist(vDB$Keywords)
	    keywords[,"ID"] <- unlist(lapply(1:nrow(vDB),
		   function(i) rep(id[i], length(vDB$Keywords[[i]]))))
	    keywords[,"Package"] <- pkg
	    hDB[[4L]] <- rbind(hDB[[4L]], keywords)
	}
    }
    hDB
}

merge.demo.index <- function(hDB, path, pkg) {
    ## Demos in the hsearch index started in R2.14.0
    if(file.exists(d_file <- file.path(path, "Meta", "demo.rds"))
       && !is.null(dDB <- readRDS(d_file))
       && nrow(dDB)) {
	## Make it look like an hDB base matrix and append it
	base <- matrix("", nrow=nrow(dDB), ncol=8)
	colnames(base) <- colnames(hDB[[1]])
	base[,"Package"] <- pkg
	base[,"LibPath"] <- path
	id <- as.character(1:nrow(dDB) + NROW(hDB[[1]]))
	base[,"ID"] <- id
	base[,"name"] <- dDB[,1]
	base[,"topic"] <- base[,"name"]
	base[,"title"] <- dDB[,2]
	base[,"Type"] <- "demo"
	hDB[[1L]] <- rbind(hDB[[1L]], base)
	aliases <- matrix("", nrow=nrow(dDB), ncol=3)
	colnames(aliases) <- colnames(hDB[[2]])
	aliases[,"Aliases"] <- base[,"name"]
	aliases[,"ID"] <- id
	aliases[,"Package"] <- pkg
	hDB[[2L]] <- rbind(hDB[[2L]], aliases)
    }
    hDB
}

## FIXME: use UTF-8, either always or optionally
## (Needs UTF-8-savvy & fast agrep, and PCRE regexps.)
help.search <-
    function(pattern, fields = c("alias", "concept", "title"),
             apropos, keyword, whatis, ignore.case = TRUE,
             package = NULL, lib.loc = NULL,
             help.db = getOption("help.db"),
             verbose = getOption("verbose"),
             rebuild = FALSE, agrep = NULL, use_UTF8 = FALSE,
             types = getOption("help.search.types")
)
{
    WINDOWS <- .Platform$OS.type == "windows"

    ### Argument handling.
    FIELDS <- c("alias", "concept", "keyword", "name", "title")
    TYPES <- c("help", "vignette", "demo")

    if (is.logical(verbose)) verbose <- 2*as.integer(verbose)
    .wrong_args <- function(args)
	gettextf("argument %s must be a single character string", sQuote(args))

    fuzzy <- agrep
    if(!missing(pattern)) {
	if(!is.character(pattern) || (length(pattern) > 1L))
	    stop(.wrong_args("pattern"), domain = NA)
	i <- pmatch(fields, FIELDS)
	if(any(is.na(i)))
	    stop("incorrect field specification")
	else
	    fields <- FIELDS[i]
    } else if(!missing(apropos)) {
	if(!is.character(apropos) || (length(apropos) > 1L))
	    stop(.wrong_args("apropos"), domain = NA)
	else {
	    pattern <- apropos
	    fields <- c("alias", "title")
	}
    } else if(!missing(keyword)) {
	if(!is.character(keyword) || (length(keyword) > 1L))
	    stop(.wrong_args("keyword"), domain = NA)
	else {
	    pattern <- keyword
	    fields <- "keyword"
	    if(is.null(fuzzy)) fuzzy <- FALSE
	}
    } else if(!missing(whatis)) {
	if(!is.character(whatis) || (length(whatis) > 1))
	    stop(.wrong_args("whatis"), domain = NA)
	else {
	    pattern <- whatis
	    fields <- "alias"
	}
    } else {
	stop("do not know what to search")
    }
    i <- pmatch(types, TYPES)
    if (any(is.na(i)))
	stop("incorrect type specification")
    else
	types <- TYPES[i]

    if(is.null(lib.loc))
	lib.loc <- .libPaths()

    if(!missing(help.db))
	warning("argument 'help.db' is deprecated")

    ### Set up the hsearch db.
    db <- eval(.hsearch_db())
    if(is.null(db))
	rebuild <- TRUE
    else if(!rebuild) {
	## Need to find out whether this has the info we need.
	## Note that when looking for packages in libraries we always
	## use the first location found.  Hence if the library search
	## path changes we might find different versions of a package.
	## Thus we need to rebuild the hsearch db in case the specified
	## library path is different from the one used when building the
	## hsearch db (stored as its "LibPaths" attribute).
	if(!identical(lib.loc, attr(db, "LibPaths")) ||
	   !all(types %in% attr(db, "Types")) ||
	   ## We also need to rebuild the hsearch db in case an existing
	   ## dir in the library path was modified more recently than
	   ## the db, as packages might have been installed or removed.
	   any(attr(db, "mtime") <
	       file.info(lib.loc[file.exists(lib.loc)])$mtime) ||
	   ## Or if the user changed the locale character type ...
	   !identical(attr(db, "ctype"), Sys.getlocale("LC_CTYPE"))
	   )
	    rebuild <- TRUE
        ## We also need to rebuild if 'packages' was used before and has
        ## changed.
        if (!is.null(package) &&
            any(! package %in% db$Base[, "Package"]))
            rebuild <- TRUE
    }
    if(rebuild) {
	if(verbose > 0L) {
            message("Rebuilding the help.search() database", " ", "...",
                    if(verbose > 1L) "...", domain = NA)
            flush.console()
        }

	if(!is.null(package)) {
	    packages_in_hsearch_db <- package
            package_paths <- NULL
	} else {
            ## local version of .packages(all.available = TRUE),
            ## recording paths
            ans <- character(0L); paths <- character(0L)
            lib.loc <- lib.loc[file.exists(lib.loc)]
            valid_package_version_regexp <-
                .standard_regexps()$valid_package_version
            for (lib in lib.loc) {
                a <- list.files(lib, all.files = FALSE, full.names = FALSE)
                for (nam in a) {
                    pfile <- file.path(lib, nam, "Meta", "package.rds")
                    if (file.exists(pfile))
                        info <- readRDS(pfile)$DESCRIPTION[c("Package", "Version")]
                    else next
                    if ( (length(info) != 2L) || any(is.na(info)) ) next
                    if (!grepl(valid_package_version_regexp, info["Version"])) next
                    ans <- c(ans, nam)
                    paths <- c(paths, file.path(lib, nam))
                }
            }
            un <- !duplicated(ans)
	    packages_in_hsearch_db <-  ans[un]
            package_paths <- paths[un]
            names(package_paths) <- ans[un]
        }

	## Create the hsearch db.
	np <- 0L
	if(verbose >= 2L) {
	    message("Packages {readRDS() sequentially}:", domain = NA)
            flush.console()
        }
        tot <- length(package_paths)
        incr <- 0L
        if(verbose && WINDOWS) {
            pb <- winProgressBar("R: creating the help.search() DB", max = tot)
            on.exit(close(pb))
        } else if(verbose == 1L) incr <- ifelse(tot > 500L, 100L, 10L)

	## Starting with R 1.8.0, prebuilt hsearch indices are available
	## in Meta/hsearch.rds, and the code to build this from the Rd
	## contents (as obtained from both new and old style Rd indices)
	## has been moved to tools:::.build_hsearch_index() which
	## creates a per-package list of base, aliases and keywords
	## information.	 When building the global index, it seems (see
	## e.g. also the code in tools:::Rdcontents()), most efficient to
	## create a list *matrix* (dbMat below), stuff the individual
	## indices into its rows, and finally create the base, alias,
	## keyword, and concept information in rbind() calls on the
	## columns.  This is *much* more efficient than building
	## incrementally.
	dbMat <- vector("list", length(packages_in_hsearch_db) * 4L)
	dim(dbMat) <- c(length(packages_in_hsearch_db), 4L)
	defunct_standard_package_names <-
	    tools:::.get_standard_package_names()$stubs

	for(p in packages_in_hsearch_db) {
            if(incr && np %% incr == 0L) {
                message(".", appendLF = FALSE, domain = NA)
                flush.console()
            }
	    np <- np + 1L
            if(verbose && WINDOWS) setWinProgressBar(pb, np)
	    if(verbose >= 2L) {
		message(" ", p, appendLF = ((np %% 5L) == 0L), domain=NA)
                flush.console()
            }
            path <- if(!is.null(package_paths)) package_paths[p]
	    else find.package(p, lib.loc, quiet = TRUE)
	    if(length(path) == 0L) {
                if(is.null(package)) next
		else stop(gettextf("could not find package %s", sQuote(p)),
                          domain = NA)
            }
	    ## Hsearch 'Meta/hsearch.rds' indices were introduced in
	    ## R 1.8.0.	 If they are missing, we really cannot use
	    ## the package (as library() will refuse to load it).
	    ## We always load hsearch.rds to establish the format,
	    ## sometimes vignettes.rds.

	    if(file.exists(hs_file <- file.path(path, "Meta", "hsearch.rds"))) {
		hDB <- readRDS(hs_file)
		if(!is.null(hDB)) {
		    ## Fill up possibly missing information.
		    if(is.na(match("Encoding", colnames(hDB[[1L]]))))
			hDB[[1L]] <- cbind(hDB[[1L]], Encoding = "")
		    nh <- NROW(hDB[[1L]])
		    hDB[[1L]] <- cbind(hDB[[1L]],
		                       Type = rep("help", nh))
		    if (nh)
		    	hDB[[1L]][, "LibPath"] <- path
		    if ("vignette" %in% types)
		    	hDB <- merge.vignette.index(hDB, path, p)
		    if ("demo" %in% types)
		    	hDB <- merge.demo.index(hDB, path, p)
		    ## Put the hsearch index for the np-th package into the
		    ## np-th row of the matrix used for aggregating.
		    dbMat[np, seq_along(hDB)] <- hDB
		} else if(verbose >= 2L) {
		    message(gettextf("package %s has empty hsearch data - strangely",
                                     sQuote(p)), domain = NA)
                    flush.console()
                }
	    }
	    else if(!is.null(package))
                warning("no hsearch.rds meta data for package ", p, domain = NA)
	}

	if(verbose >= 2L)  {
	    message(ifelse(np %% 5L == 0L, "\n", "\n\n"),
                    sprintf("Built dbMat[%d,%d]", nrow(dbMat), ncol(dbMat)),
                    domain = NA)
            flush.console()
            ## DEBUG save(dbMat, file="~/R/hsearch_dbMat.rda", compress=TRUE)
        }

	## workaround methods:::rbind() misbehavior:
	if(.isMethodsDispatchOn()) {
	    bind_was_on <- methods:::bind_activation(FALSE)
	    if(bind_was_on) on.exit(methods:::bind_activation(TRUE))
	}

	## Create the global base, aliases, keywords and concepts tables
	## via calls to rbind() on the columns of the matrix used for
	## aggregating.
	db <- list(Base     = do.call("rbind", dbMat[, 1]),
		   Aliases  = do.call("rbind", dbMat[, 2]),
		   Keywords = do.call("rbind", dbMat[, 3]),
		   Concepts = do.call("rbind", dbMat[, 4]))
	if(is.null(db$Concepts))
	    db$Concepts <-
		matrix(character(), ncol = 3L,
		       dimnames = list(NULL,
		       c("Concepts", "ID", "Package")))
	## Make the IDs globally unique by prefixing them with the
	## number of the package in the global index.
	for(i in which(sapply(db, NROW) > 0L)) {
	    db[[i]][, "ID"] <-
		paste(rep.int(seq_along(packages_in_hsearch_db),
			      sapply(dbMat[, i], NROW)),
		      db[[i]][, "ID"],
		      sep = "/")
	}
	## And maybe re-encode ...
	if(!identical(Sys.getlocale("LC_CTYPE"), "C")) {
	    if(verbose >= 2L) {
                message("reencoding ...", appendLF=FALSE, domain = NA)
                flush.console()
            }
	    encoding <- db$Base[, "Encoding"]
            target <- ifelse(use_UTF8 && !l10n_info()$`UTF-8`, "UTF-8", "")
	    ## As iconv is not vectorized in the 'from' argument, loop
	    ## over groups of identical encodings.
	    for(enc in unique(encoding)) {
                if(enc != target) next
		IDs <- db$Base[encoding == enc, "ID"]
		for(i in seq_along(db)) {
		    ind <- db[[i]][, "ID"] %in% IDs
		    db[[i]][ind, ] <- iconv(db[[i]][ind, ], enc, "")
		}
	    }
	    if(verbose >= 2L) {
                message(" ", "done", domain = NA)
                flush.console()
            }
	}
	bad_IDs <-
	    unlist(sapply(db,
			  function(u)
			  u[rowSums(is.na(nchar(u, "c", TRUE))) > 0, "ID"]))
        ## FIXME: drop this fallback
	if(length(bad_IDs)) { ## try latin1
            for(i in seq_along(db)) {
                ind <- db[[i]][, "ID"] %in% bad_IDs
                db[[i]][ind, ] <- iconv(db[[i]][ind, ], "latin1", "")
            }
            bad_IDs <-
                unlist(sapply(db,
                              function(u)
                              u[rowSums(is.na(nchar(u, "c", TRUE))) > 0, "ID"]))
        }
	## If there are any invalid multi-byte character data
	## left, we simple remove all Rd objects with at least one
	## invalid entry, and warn.
        if(length(bad_IDs)) {
	    warning("removing all entries with invalid multi-byte character data")
	    for(i in seq_along(db)) {
		ind <- db[[i]][, "ID"] %in% bad_IDs
		db[[i]] <- db[[i]][!ind, ]
	    }
	}

        if(verbose >= 2L) {
            message("saving the database ...", appendLF=FALSE, domain = NA)
            flush.console()
        }
        attr(db, "LibPaths") <- lib.loc
        attr(db, "mtime") <- Sys.time()
        attr(db, "ctype") <- Sys.getlocale("LC_CTYPE")
        attr(db, "Types") <- unique(c("help", types))
        .hsearch_db(db)
        if(verbose >= 2L) {
            message(" ", "done", domain = NA)
            flush.console()
        }
        if(verbose > 0L) {
            message("... database rebuilt", domain = NA)
            if(WINDOWS) {
                close(pb)
                on.exit() # clear closing of progress bar
            }
            flush.console()
        }
    }

    ### Matching.
    if(verbose >= 2L) {
	message("Database of ",
                NROW(db$Base), " help objects (",
                NROW(db$Aliases), " aliases, ",
                NROW(db$Concepts), " concepts, ",
                NROW(db$Keywords), " keywords)",
                domain = NA)
        flush.console()
    }
    if(!is.null(package)) {
	## Argument 'package' was given.  Need to check that all given
	## packages exist in the db, and only search the given ones.
	pos_in_hsearch_db <-
	    match(package, unique(db$Base[, "Package"]), nomatch = 0L)
        ## This should not happen for R >= 2.4.0
	if(any(pos_in_hsearch_db) == 0L)
	    stop(gettextf("no information in the database for package %s: need 'rebuild = TRUE'?",
			  sQuote(package[pos_in_hsearch_db == 0][1L])),
                 domain = NA)
	db <-
	    lapply(db,
		   function(x) {
		       x[x[, "Package"] %in% package, , drop = FALSE]
		   })
    }

    ## Subset to the requested help types
    db$Base <- db$Base[db$Base[,"Type"] %in% types,,drop=FALSE]

    ## <FIXME>
    ## No need continuing if there are no objects in the data base.
    ## But shouldn't we return something of class "hsearch"?
    if(!length(db$Base)) return(invisible())
    ## </FIXME>

    ## If agrep is NULL (default), we want to use fuzzy matching iff
    ## 'pattern' contains no characters special to regular expressions.
    ## We use the following crude approximation: if pattern contains
    ## only alphanumeric characters or whitespace or a '-', it is taken
    ## 'as is', and fuzzy matching is used unless turned off explicitly,
    ## or pattern has very few (currently, less than 5) characters.
    if(is.null(fuzzy) || is.na(fuzzy))
	fuzzy <-
	    (grepl("^([[:alnum:]]|[[:space:]]|-)+$", pattern)
	     && (nchar(pattern, type="c") > 4L))
    if(is.logical(fuzzy)) {
	if(fuzzy)
	    max.distance <- 0.1
    }
    else if(is.numeric(fuzzy) || is.list(fuzzy)) {
	max.distance <- fuzzy
	fuzzy <- TRUE
    }
    else
	stop("incorrect 'agrep' specification")

    searchFun <- function(x) {
	if(fuzzy)
	    agrep(pattern, x, ignore.case = ignore.case,
		  max.distance = max.distance)
	else
	    grep(pattern, x, ignore.case = ignore.case, perl = use_UTF8)
    }
    dbBase <- db$Base
    searchDbField <- function(field) {
	switch(field,
	       alias = {
		   aliases <- db$Aliases
		   match(aliases[searchFun(aliases[, "Aliases"]),
				 "ID"],
			 dbBase[, "ID"])
	       },
	       concept = {
		   concepts <- db$Concepts
		   match(concepts[searchFun(concepts[, "Concepts"]),
				  "ID"],
			 dbBase[, "ID"])
	       },

	       keyword = {
		   keywords <- db$Keywords
		   match(keywords[searchFun(keywords[, "Keywords"]),
				  "ID"],
			 dbBase[, "ID"])
	       },
	       searchFun(db$Base[, field]))
    }

    i <- NULL
    for(f in fields) i <- c(i, searchDbField(f))
    db <- dbBase[sort(unique(i)),
		 c("topic", "title", "Package", "LibPath", "name", "Type"),
		 drop = FALSE]
    if(verbose>= 2L) {
        message(sprintf(ngettext(NROW(db),
                                 "matched %d object.",
                                 "matched %d objects."),
                        NROW(db)),
                domain = NA)
        flush.console()
    }

    ## Retval.
    y <- list(pattern = pattern, fields = fields,
	      type = if(fuzzy) "fuzzy" else "regexp",
	      agrep = agrep,
	      ignore.case = ignore.case, types = types,
	      package = package, lib.loc = lib.loc,
	      matches = db)
    class(y) <- "hsearch"
    y
}

## this extra indirection allows the Mac GUI to replace this
## yet call the printhsearchInternal function.
print.hsearch <- function(x, ...)
    printhsearchInternal(x, ...)


printhsearchInternal  <- function(x, ...)
{
    help_type <- getOption("help_type", default="text")
    types <- x$types
    if (help_type == "html") {
        browser <- getOption("browser")
	if (tools:::httpdPort == 0L) tools::startDynamicHelp()
	if (tools:::httpdPort > 0L) {
	    url <- paste0("http://127.0.0.1:", tools:::httpdPort,
                      "/doc/html/Search?pattern=", tools:::escapeAmpersand(x$pattern),
                      # Only encode non-default values
                      if (!("title" %in% x$fields)) "&title=0",
                      if ("keyword" %in% x$fields) "&keyword=1",
                      if (!("alias" %in% x$fields)) "&alias=0",
                      if (!("concept" %in% x$fields)) "&concept=0",
                      if ("name" %in% x$fields) "&name=1",
                      if (!is.null(x$agrep)) paste0("&agrep=", x$agrep),
                      if (!x$ignore.case) "&ignore.case=0",
                      if (!identical(types, getOption("help.search.types")))
			 paste0("&types=", paste(types, collapse=";")),
                      if (!is.null(x$package))
			 paste0("&package=", paste(x$package, collapse=";")),
                      if (!identical(x$lib.loc, .libPaths()))
			 paste0("&lib.loc=", paste(x$lib.loc, collapse=";")))
            browseURL(url, browser)
            return(invisible(x))
        }
    }
    hfields <- paste(x$fields, collapse = " or ")
    vfieldnames <- c(alias = "name", concept="keyword", keyword=NA,
                     name="name", title="title")
    vfieldnames <- vfieldnames[x$fields]
    vfields <- paste(unique(vfieldnames[!is.na(vfieldnames)]), collapse = " or ")
    dfieldnames <- c(alias = "name", concept=NA, keyword=NA,
                     name = "name", title = "title")
    dfieldnames <- dfieldnames[x$fields]
    dfields <- paste(unique(dfieldnames[!is.na(dfieldnames)]), collapse = " or ")
    fields <- list(help=hfields, vignette=vfields, demo=dfields)
    matchtype <- switch(x$type, fuzzy = "fuzzy", "regular expression")
    typenames <- c(vignette = "Vignettes", help = "Help files", demo="Demos")
    db <- x$matches
    if(NROW(db) == 0) {
    	typenames <- paste(tolower(typenames[types]), collapse=" or ")
	writeLines(strwrap(paste("No", typenames,  "found with", fields$help,
				 "matching", sQuote(x$pattern),
				 "using", matchtype, "matching.")))
        return(invisible(x))
    }

    outFile <- tempfile()
    outConn <- file(outFile, open = "w")
    typeinstruct <- c(vignette = paste("Type 'vignette(\"FOO\", package=\"PKG\")' to",
				       "inspect entries 'PKG::FOO'."),
                      help = paste("Type '?PKG::FOO' to",
				       "inspect entries 'PKG::FOO',",
				       "or 'TYPE?PKG::FOO' for entries like",
				       "'PKG::FOO-TYPE'."),
		      demo = paste("Type 'demo(PKG::FOO)' to",
				       "run demonstration 'PKG::FOO'."))

    for (type in types) {
	if(NROW(dbtemp <- db[db[,"Type"] == type,,drop=FALSE]) > 0) {
	    writeLines(c(strwrap(paste(typenames[type], "with", fields[[type]],
				       "matching", sQuote(x$pattern),
				       "using", matchtype, "matching:")),
			 "\n"),
		       outConn)
	    dbnam <- paste0(dbtemp[, "Package"], "::", dbtemp[ , "topic"])
	    dbtit <- paste0(dbtemp[ , "title"])
	    writeLines(formatDL(dbnam, dbtit), outConn)
	    writeLines(c("\n",
			 strwrap(typeinstruct[type]),
			 "\n\n"),
		       outConn)
	}
    }
    close(outConn)
    file.show(outFile, delete.file = TRUE)
    invisible(x)
}
#  File src/library/utils/R/help.start.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

help.start <-
    function (update = FALSE, gui = "irrelevant",
              browser = getOption("browser"), remote = NULL)
{
    WINDOWS <- .Platform$OS.type == "windows"
    if (!WINDOWS) {
        ## should always be set, but might be empty
        if (!is.function(browser) &&
            (length(browser) != 1L || !is.character(browser) || !nzchar(browser)))
            stop("invalid browser name, check options(\"browser\").")
    }
    home <- if (is.null(remote)) {
        if (tools:::httpdPort == 0L) tools::startDynamicHelp()
        if (tools:::httpdPort > 0L) {
            if (update) make.packages.html(temp = TRUE)
            paste0("http://127.0.0.1:", tools:::httpdPort)
        } else stop("help.start() requires the HTTP server to be running",
                    call. = FALSE)
    } else remote
    url <- paste0(home, "/doc/html/index.html")

    ## FIXME: maybe these should use message()?
    if (WINDOWS) {
        cat(gettextf("If nothing happens, you should open\n'%s' yourself\n", url))
    } else if (is.character(browser)) {
        writeLines(strwrap(gettextf("If the browser launched by '%s' is already running, it is *not* restarted, and you must switch to its window.",
                                    browser),
                           exdent = 4L))
        writeLines(gettext("Otherwise, be patient ..."))
    }
    browseURL(url, browser = browser)
    invisible()
}

browseURL <- function(url, browser = getOption("browser"), encodeIfNeeded=FALSE)
{
    WINDOWS <- .Platform$OS.type == "windows"

    if (!is.character(url) || length(url) != 1L|| !nzchar(url))
        stop("'url' must be a non-empty character string")
    if(identical(browser, "false")) return(invisible())
    if(WINDOWS && is.null(browser)) return(shell.exec(url))
    else if (is.function(browser))
        return(invisible(browser(if(encodeIfNeeded) URLencode(url) else url)))

   if (!is.character(browser) || length(browser) != 1L || !nzchar(browser))
        stop("'browser' must be a non-empty character string")
    if (WINDOWS) {
        ## No shell used, but spaces are possible
        return(system(paste0('"', browser, '" ',
                             if(encodeIfNeeded) URLencode(url) else url),
                      wait = FALSE))
    }

    ## Unix-alike, character "browser"
    if (.Platform$GUI == "AQUA" ||
        length(grep("^(localhost|):", Sys.getenv("DISPLAY"))) )
      isLocal <- TRUE
    else
      isLocal <- FALSE

    ## escape characters.  ' can occur in URLs, so we must use " to
    ## delimit the URL.  We need to escape $, but "`\ do not occur in
    ## valid URLs (RFC 2396, on the W3C site).
    .shQuote <- function(string)
        paste0('"', gsub("\\$", "\\\\$", string), '"')
    quotedUrl <- .shQuote(if(encodeIfNeeded) URLencode(url) else url)
    remoteCmd <- if (isLocal)
        switch(basename(browser),
               "gnome-moz-remote" =, "open" = quotedUrl,
               "galeon" = paste("-x", quotedUrl),
               "kfmclient" = paste("openURL", quotedUrl),
               "mozilla" =, "opera" =, "firefox" = {
                   paste0("-remote \"openURL(",
                         ## Quote ',' and ')' ...
                         gsub("([,)$])", "%\\1", url), ")\"")
               }, quotedUrl)
    else quotedUrl
    system(paste(browser, remoteCmd, "> /dev/null 2>&1 ||",
                 browser, quotedUrl, "&"))
}
#  File src/library/utils/R/history.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

loadhistory <- function(file = ".Rhistory")
    invisible(.External2(C_loadhistory, file))

savehistory <- function(file = ".Rhistory")
    invisible(.External2(C_savehistory, file))

history <- function(max.show = 25, reverse = FALSE, pattern, ...)
{
    file1 <- tempfile("Rrawhist")
    savehistory(file1)
    rawhist <- readLines(file1)
    unlink(file1)
    if(!missing(pattern))
        rawhist <- unique(grep(pattern, rawhist, value = TRUE, ...))
    nlines <- length(rawhist)
    if(nlines) {
        inds <- max(1, nlines-max.show):nlines
        if(reverse) inds <- rev(inds)
    } else inds <- integer()
    file2 <- tempfile("hist")
    writeLines(rawhist[inds], file2)
    file.show(file2, title = "R History", delete.file = TRUE)
}

timestamp <- function(stamp = date(), prefix = "##------ ",
                      suffix = " ------##", quiet = FALSE)
{
    stamp <- paste0(prefix, stamp, suffix)
    .External2(C_addhistory, stamp)
    if (!quiet) cat(stamp, sep = "\n")
    invisible(stamp)
}
#  File src/library/utils/R/iconv.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


## If you were wondering what these language codes stand for, see
## ftp://ftp.ilog.fr/pub/Users/haible/utf8/ISO_639
localeToCharset <- function(locale = Sys.getlocale("LC_CTYPE"))
{
    guess <- function(en)
    {
        if(en %in% c("aa", "af", "an", "br", "ca", "da", "de", "en",
                         "es", "et", "eu", "fi", "fo", "fr", "ga", "gl",
                         "gv", "id", "is", "it", "kl", "kw", "ml", "ms",
                         "nb", "nn", "no", "oc", "om", "pt", "so", "sq",
                         "st", "sv", "tl", "uz", "wa", "xh", "zu"))
                return("ISO8859-1")
        if(en %in% c("bs", "cs", "hr", "hu", "pl", "ro", "sk", "sl"))
            return("ISO8859-2")
        if(en %in% "mt") return("ISO8859-3")
            if(en %in% c("mk", "ru")) return("ISO8859-5")
        if(en %in% "ar") return("ISO8859-6")
        if(en %in% "el") return("ISO8859-7")
        if(en %in% c("he", "iw", "yi")) return("ISO8859-8")
        if(en %in% "tr") return("ISO8859-9")
        if(en %in% "lg") return("ISO8859-10")
        if(en %in% c("lt", "lv", "mi")) return("ISO8859-13")
        if(en %in% "cy") return("ISO8859-14")
        if(en %in% "uk") return("KOI8-U")
        if(en %in% "ja") return("EUC-JP")
        if(en %in% "ko") return("EUC-KR")
        if(en %in% "th") return("TIS-620")
        if(en %in% "tg") return("KOI8-T")
        if(en %in% "ka") return("GEORGIAN-PS")
        if(en %in% "kk") return("PT154")
        ## not safe to guess for zh
        return(NA_character_)
    }
    if(locale %in% c("C", "POSIX")) return("ASCII")
    if(.Platform$OS.type == "windows") {
        x <- strsplit(locale, ".", fixed=TRUE)[[1L]]
        if(length(x) != 2) return(NA_character_)
        ## PUTTY suggests mapping Windows code pages as
        ## 1250 -> ISO 8859-2
        ## 1251 -> KOI8-U
        ## 1252 -> ISO 8859-1
        ## 1253 -> ISO 8859-7
        ## 1254 -> ISO 8859-9
        ## 1255 -> ISO 8859-8
        ## 1256 -> ISO 8859-6
        ## 1257 -> ISO 8859-13
        switch(x[2L],
              # this is quite wrong "1250" = return("ISO8859-2"),
              # this is quite wrong "1251" = return("KOI8-U"),
               "1252" = return("ISO8859-1"),
              # "1253" = return("ISO8859-7"),
              # "1254" = return("ISO8859-9"),
              # "1255" = return("ISO8859-8"),
              # "1256" = return("ISO8859-6"),
               "1257" = return("ISO8859-13")
               )
        return(paste0("CP", x[2L]))
    } else {
        ## Assume locales are like  en_US[.utf8[@euro]]
        x <- strsplit(locale, ".", fixed=TRUE)[[1L]]
        enc <- if(length(x) == 2) gsub("@.*$o", "", x[2L]) else ""
        # AIX uses UTF-8, OS X utf-8
        if(toupper(enc) == "UTF-8") enc <- "utf8"
        if(nzchar(enc) && enc != "utf8") {
            enc <- tolower(enc)
            known <-
                c("ISO8859-1", "ISO8859-2", "ISO8859-3", "ISO8859-6",
                  "ISO8859-7", "ISO8859-8", "ISO8859-9", "ISO8859-10",
                  "ISO8859-13", "ISO8859-14", "ISO8859-15",
                  "CP1251", "CP1255", "EUC-JP", "EUC-KR", "EUC-TW",
                  "GEORGIAN-PS", "KOI8-R",  "KOI8-U", "TCVN",
                  "BIG5" , "GB2312", "GB18030", "GBK",
                  "TIS-620", "SHIFT_JIS", "GB2312", "BIG5-HKSCS")
            names(known) <-
                c("iso88591", "iso88592", "iso88593", "iso88596",
                  "iso88597", "iso88598", "iso88599", "iso885910",
                  "iso885913", "iso885914", "iso885915",
                  "cp1251", "cp1255", "eucjp", "euckr", "euctw",
                  "georgianps", "koi8r", "koi8u", "tcvn",
                  "big5" , "gb2312", "gb18030", "gbk",
                  "tis-620", "sjis", "eucn", "big5-hkscs")
	    if (length(grep("darwin",R.version$os))) {
	        k <- c(known, "ISO8859-1", "ISO8859-2", "ISO8859-4",
		  "ISO8859-7", "ISO8859-9", "ISO8859-13", "ISO8859-15",
		  "KOI8-U", "KOI8-R", "PT154", "ASCII", "ARMSCII-8",
		  "ISCII-DEV", "BIG5-HKCSC")
		names(k) <- c(names(known), "iso8859-1", "iso8859-2", "iso8859-4",
		  "iso8859-7", "iso8859-9", "iso8859-13", "iso8859-15",
		  "koi8-u", "koi8-r", "pt154", "us-ascii", "armscii-8",
		  "iscii-dev", "big5hkscs")
		known <- k
            }
	    if(enc %in% names(known)) return(unname(known[enc]))
            if(length(grep("^cp-", enc)))  # old Linux
                return(sub("cp-([0-9]+)", "CP\\1", enc))
            if(enc == "EUC") {
                ## let's hope it is a ll_* name.
                if(length(grep("^[[:alpha:]]{2}_", x[1L], perl = TRUE))) {
                    ll <- substr(x[1L], 1L, 2L)
                    return(switch(ll, "jp"="EUC-JP", "kr"="EUC-KR",
                                  "zh"="GB2312"))
                }
            }
        }
	## on Darwin all real locales w/o encoding are UTF-8
	## HOWEVER! unlike the C code, we cannot filter out
	## invalid locales, so it will be wrong for non-supported
	## locales (why is this duplicated in R code anyway?)
	if (length(grep("darwin", R.version$os))) return("UTF-8")
        ## let's hope it is a ll_* name.
        if(length(grep("^[[:alpha:]]{2}_", x[1L], perl = TRUE))) {
            ll <- substr(x[1L], 1L, 2L)
            if(enc == "utf8") return(c("UTF-8", guess(ll)))
            else return(guess(ll))
        }
        return(NA_character_)
    }
}
#  File src/library/utils/R/indices.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

packageDescription <- function(pkg, lib.loc=NULL, fields=NULL, drop=TRUE,
			       encoding = "")
{
    retval <- list()
    if(!is.null(fields)){
        fields <- as.character(fields)
        retval[fields] <- NA
    }

    ## If the NULL default for lib.loc is used, the loaded packages are
    ## searched before the libraries.
    pkgpath <-
	if(is.null(lib.loc)) {
	    if(pkg == "base")
		file.path(.Library, "base")
	    else if((envname <- paste0("package:", pkg))
		    %in% search()) {
		pp <- attr(as.environment(envname), "path")
		## could be NULL if a perverse user has been naming
		## environments to look like packages.
	    } else if(pkg %in% loadedNamespaces())
		## correct path for a loaded (not attached) namespace:
		getNamespaceInfo(pkg, "path")
	}
    if(is.null(pkgpath)) pkgpath <- ""

    if(pkgpath == "") {
        libs <- if(is.null(lib.loc)) .libPaths() else lib.loc
        for(lib in libs)
            if(file.access(file.path(lib, pkg), 5) == 0L) {
                pkgpath <- file.path(lib, pkg)
                break
            }
    }

    if(pkgpath == "") {
        warning(gettextf("no package '%s' was found", pkg), domain = NA)
        return(NA)
    }

    ## New in 2.7.0: look for installed metadata first.
    ## We always need to be able to drop back to the file as this
    ## is used during package installation.

    if(file.exists(file <- file.path(pkgpath, "Meta", "package.rds"))) {
        desc <- readRDS(file)$DESCRIPTION
        if(length(desc) < 1)
            stop(gettextf("metadata of package '%s' is corrupt", pkg),
                 domain = NA)
        desc <- as.list(desc)
    } else if(file.exists(file <- file.path(pkgpath,"DESCRIPTION"))) {
        dcf <- read.dcf(file=file)
        if(NROW(dcf) < 1L)
            stop(gettextf("DESCRIPTION file of package '%s' is corrupt", pkg),
                 domain = NA)
        desc <- as.list(dcf[1,])
    } else file <- ""

    if(file != "") {
        ## read the Encoding field if any
        enc <- desc[["Encoding"]]
        if(!is.null(enc) && !is.na(encoding)) {
            ## Determine encoding and re-encode if necessary and possible.
            if (missing(encoding) && Sys.getlocale("LC_CTYPE") == "C")
                encoding <- "ASCII//TRANSLIT"
            ## might have an invalid encoding ...
            newdesc <- try(lapply(desc, iconv, from = enc, to = encoding))
            if(!inherits(newdesc, "try-error")) desc <- newdesc
            else
                warning("'DESCRIPTION' file has an 'Encoding' field and re-encoding is not possible", call. = FALSE)
        }
        if(!is.null(fields)){
            ok <- names(desc) %in% fields
            retval[names(desc)[ok]] <- desc[ok]
        }
        else
            retval[names(desc)] <- desc
    }

    if((file == "") || (length(retval) == 0)){
        warning(gettextf("DESCRIPTION file of package '%s' is missing or broken", pkg), domain = NA)
        return(NA)
    }

    if(drop & length(fields) == 1L)
        return(retval[[1L]])

    class(retval) <- "packageDescription"
    if(!is.null(fields)) attr(retval, "fields") <- fields
    attr(retval, "file") <- file
    retval
}


print.packageDescription <-
    function(x, abbrCollate = 0.8 * getOption("width"), ...)
{
    xx <- x
    xx[] <- lapply(xx, function(x) if(is.na(x)) "NA" else x)
    if(abbrCollate > 0 && any(names(xx) == "Collate")) {
        ## trim a long "Collate" field -- respecting word boundaries
	wrds <- strsplit(xx$Collate,"[ \n]")[[1L]]
	k <- which.max(cumsum(nchar(wrds)) > abbrCollate) - 1L
	xx$Collate <- paste(c(wrds[seq_len(k)], "....."), collapse=" ")
    }
    write.dcf(as.data.frame.list(xx, optional = TRUE))
    cat("\n-- File:", attr(x, "file"), "\n")
    if(!is.null(attr(x, "fields"))){
        cat("-- Fields read: ")
        cat(attr(x, "fields"), sep = ", ")
        cat("\n")
    }
    invisible(x)
}

# Simple convenience functions

maintainer <- function(pkg)
{
    force(pkg)
    desc <- packageDescription(pkg)
    if(is.list(desc)) gsub("\n", " ", desc$Maintainer, fixed = TRUE)
    else NA_character_
}

packageVersion <- function(pkg, lib.loc = NULL)
{
    res <- suppressWarnings(packageDescription(pkg, lib.loc=lib.loc,
                                               fields = "Version"))
    if (!is.na(res)) package_version(res) else
    stop(gettextf("package %s not found", sQuote(pkg)), domain = NA)
}

## used with firstOnly = TRUE for example()
## used with firstOnly = FALSE in help()
index.search <- function(topic, paths, firstOnly = FALSE)
{
    res <- character()
    for (p in paths) {
        if(file.exists(f <- file.path(p, "help", "aliases.rds")))
            al <- readRDS(f)
        else if(file.exists(f <- file.path(p, "help", "AnIndex"))) {
            ## aliases.rds was introduced before 2.10.0, as can phase this out
            foo <- scan(f, what = list(a="", b=""), sep = "\t", quote = "",
                        na.strings = "", quiet = TRUE)
            al <- structure(foo$b, names = foo$a)
        } else next
        f <- al[topic]
        if(is.na(f)) next
        res <- c(res, file.path(p, "help", f))
        if(firstOnly) break
    }
    res
}

print.packageIQR <-
function(x, ...)
{
    db <- x$results
    ## Split according to Package.
    out <- if(nrow(db) == 0L)
         NULL
    else
        lapply(split(1 : nrow(db), db[, "Package"]),
               function(ind) db[ind, c("Item", "Title"),
                                drop = FALSE])
    outFile <- tempfile("RpackageIQR")
    outConn <- file(outFile, open = "w")
    first <- TRUE
    for(pkg in names(out)) {
        writeLines(paste0(ifelse(first, "", "\n"), x$title,
                          " in package ", sQuote(pkg), ":\n"),
                   outConn)
        writeLines(formatDL(out[[pkg]][, "Item"],
                            out[[pkg]][, "Title"]),
                   outConn)
        first <- FALSE
    }
    if(first) {
        close(outConn)
        unlink(outFile)
        writeLines(paste("no", tolower(x$title), "found"))
        if(!is.null(x$footer))
            writeLines(c("", x$footer))
    }
    else {
        if(!is.null(x$footer))
            writeLines(c("\n", x$footer), outConn)
        close(outConn)
        file.show(outFile, delete.file = TRUE,
                  title = paste("R", tolower(x$title)))
    }
    invisible(x)
}
#  File src/library/utils/R/linkhtml.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


make.packages.html <-
    function(lib.loc = .libPaths(), temp = FALSE, verbose = TRUE,
             docdir = R.home("doc"))
{
    add_lib_index <- function(libs)
    {
        cat('<div align="left">\n<ul>\n', file = out)
        for (i in seq_along(libs)) {
            nm <- libs[i]
            if (nm == .Library) {
                cat('<li>Contents of the <a href="#lib-', i, '">',
                    'standard</a> library</li>\n', sep = "", file = out)
            } else {
                cat('<li>Contents of <a href="#lib-', i, '">', nm,
                    '</a></li>\n', sep = "", file = out)
            }
        }
        cat("</ul>\n</div>\n", file = out)
    }

    WINDOWS <- .Platform$OS.type == "windows"
    f.tg <- if (temp) {
        dir.create(file.path(tempdir(), ".R/doc/html"), recursive = TRUE,
                   showWarnings = FALSE)
        file.path(tempdir(), ".R/doc/html/packages.html")
    } else file.path(docdir, "html", "packages.html")
    op <- file.path(tempdir(), ".R/doc/html/libPaths.rds")
    if (temp && file.exists(f.tg) && file.exists(op)) {
        ## check if we can avoid remaking it.
        if(identical(lib.loc, readRDS(op))) {
            dates <- file.info(c(f.tg, lib.loc))$mtime
            if(which.max(dates) == 1L) return(TRUE)
        }
    }
    if (!file.create(f.tg)) {
        warning("cannot update HTML package index")
        return(FALSE)
    }
    if (verbose) {
        message("Making 'packages.html' ...", appendLF = FALSE, domain = NA)
        flush.console()
    }
    file.append(f.tg,
                file.path(R.home("doc"), "html", "packages-head-utf8.html"))
    out <- file(f.tg, open = "a")
    on.exit(close(out))
    if(WINDOWS) {
        rh <- chartr("\\", "/", R.home())
        drive <- substring(rh, 1L, 2L)
    }
    ## find out how many
    pkgs <- vector("list", length(lib.loc))
    names(pkgs) <- lib.loc
    for (lib in lib.loc) {
        pg <- .packages(all.available = TRUE, lib.loc = lib)
        pkgs[[lib]] <- pg[order(toupper(pg), pg)]
    }
    if (WINDOWS) {
        tot <- sum(sapply(pkgs, length))
        if(verbose) {
            pb <- winProgressBar("R: creating packages.html", max = tot)
            on.exit(close(pb), add = TRUE)
        }
        npkgs <- 0L
    }
    ## If there is more than one lib, have an index at the top and bottom
    if (length(lib.loc) > 1L) add_lib_index(lib.loc)
    for (ii in seq_along(lib.loc)) {
        lib <- lib.loc[ii]
        libname <-
            if (identical(lib, .Library)) "the standard library" else if (WINDOWS) chartr("/", "\\", lib) else lib
        cat("<p><h3 id=\"lib-",ii,"\">Packages in ", libname, "</h3>\n", sep = "", file = out)
        lib0 <- "../../library"
        if (!temp) {
            if (WINDOWS) {
                ## use relative indexing for .Library
                ## perhaps other site libraries
                if (is.na(pmatch(rh, lib))) {
                    lib0 <- if(substring(lib, 2L, 2L) != ":")
                        paste0(drive, lib) else lib
                    lib0 <- paste0("file:///", URLencode(lib0))
                }
            } else {
                if (lib != .Library)
                    lib0 <- paste0("file:///", URLencode(lib))
            }
        }
        pg <- pkgs[[lib]]
        use_alpha <- (length(pg) > 100)
        first <- toupper(substr(pg, 1, 1))
        nm <- sort(names(table(first)))
        if(use_alpha) {
            writeLines("<p align=\"center\">", out)
            writeLines(paste0("<a href=\"#pkgs-", nm, "\">", nm, "</a>"), out)
            writeLines("</p>\n", out)
        }
        cat('<p><table width="100%" summary="R Package list>\n', file = out)
        for (a in nm) {
            if(use_alpha)
                cat("<tr id=\"pkgs-", a, "\"/>\n", sep = "", file = out)
            for (i in pg[first == a]) {
                title <- packageDescription(i, lib.loc = lib, fields = "Title",
                                            encoding = "UTF-8")
                if (is.na(title)) title <- "-- Title is missing --"
                cat('<tr align="left" valign="top" id="lib-"', i, '">\n',
                    '<td width="25%"><a href="', lib0, '/', i,
                    '/html/00Index.html">', i, "</a></td><td>", title,
                    "</td></tr>\n", file = out, sep = "")
                if (WINDOWS) {
                    npkgs <- npkgs + 1L
                    if(verbose) setWinProgressBar(pb, npkgs)
                }
            }
        }
        cat("</table>\n\n", file=out)
    }
    if (length(lib.loc) > 1L) add_lib_index(lib.loc)
    cat("</body></html>\n", file=out)
    if (verbose) { message(" ", "done"); flush.console() }
    if (temp) saveRDS(lib.loc, op)
    invisible(TRUE)
}
#  File src/library/utils/R/menu.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

menu <- function(choices, graphics = FALSE, title = NULL)
{
    if(!interactive()) stop("menu() cannot be used non-interactively")
    if(isTRUE(graphics)) {
        if(.Platform$OS.type == "windows" || .Platform$GUI == "AQUA"
           ## Tk might not require X11 on Mac OS X, but if DISPLAY is set
           ## this will work for Aqua Tcl/Tk.
           ## OTOH, we do want to check Tk works!
           || (capabilities("tcltk") && capabilities("X11") &&
               suppressWarnings(tcltk:::.TkUp))) {
            res <- select.list(choices, multiple = FALSE, title = title,
                               graphics = TRUE)
            return(match(res, choices, nomatch = 0L))
        }
    }
    nc <- length(choices)
    if(length(title) && nzchar(title[1L])) cat(title[1L], "\n")
    op <- paste0(format(seq_len(nc)), ": ", choices)
    if(nc > 10L) {
        fop <- format(op)
        nw <- nchar(fop[1L], "w") + 2
        ncol <- getOption("width") %/% nw  # might be 0
        if(ncol > 1L)
            op <- paste0(fop, c(rep("  ", ncol - 1), "\n"), collapse="")
        cat("", op, "", sep="\n")
    } else cat("", op, "", sep="\n")
    repeat {
	ind <- .Call(C_menu, as.character(choices))
	if(ind <= nc) return(ind)
	cat(gettext("Enter an item from the menu, or 0 to exit\n"))
    }
}
#  File src/library/utils/R/mirrorAdmin.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

mirror2html <- function(mirrors = NULL, file="mirrors.html",
                        head = "mirrors-head.html",
                        foot = "mirrors-foot.html")
{
    if(is.null(mirrors)){
        mirrors <- getCRANmirrors(all=FALSE, local.only=TRUE)
    }
    mirrors$Host <- gsub("&", "&amp;", mirrors$Host)
    z <- NULL
    if(file.exists(head)) z <- readLines(head)
    z <- c(z, "<dl>")
    for(country in unique(mirrors$Country)) {
        m <- mirrors[mirrors$Country == country,]
        z <- c(z, paste0("<dt>", country, "</dt>"),
               "<dd>",
               sprintf("<table border=0 width=90%% summary=\"%s\">",
                       country))
        for(k in seq_len(nrow(m))) {
            z <- c(z, "<tr>",
                   "<td width=45%>",
                   sprintf("<a href=\"%s\" target=\"_top\">%s</a>",
                           m$URL[k], m$URL[k]),
                   "</td>\n",
                   "<td>", m$Host[k], "</td>",
                   "</tr>")
        }
        z <- c(z, "</table>", "</dd>")
    }
    z <- c(z, "</dl>")
    if(file.exists(foot)) z <- c(z, readLines(foot))
    if(file!="") writeLines(z, file)
    invisible(z)
}

checkCRAN <- function(method)
{
    master <- available.packages(contrib.url("http://cran.R-project.org"),
                                 method=method)
    m <- getCRANmirrors()
    z <- list()
    for(url in as.character(m$URL))
        z[[url]] = available.packages(contrib.url(url), method=method)
    lapply(z, function(a) all.equal(a, master))
}






#  File src/library/utils/R/modifyList.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

### Originates from Deepayan Sarkar as  updateList() from 'lattice' package

modifyList <- function(x, val, keep.null = FALSE)
{
    stopifnot(is.list(x), is.list(val))
    xnames <- names(x)
    vnames <- names(val)
    ## Will not update unnamed components.  FIXME: What if names are repeated? Warn?
    vnames <- vnames[vnames != ""]
    if (keep.null) {
        for (v in vnames) {
            x[v] <-
                if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]]))
                    list(modifyList(x[[v]], val[[v]], keep.null = keep.null))
                else val[v]
        }
    }
    else { 
        for (v in vnames) {
            x[[v]] <-
                if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]]))
                    modifyList(x[[v]], val[[v]], keep.null = keep.null)
                else val[[v]]
        }
    }
    x
}
#  File src/library/utils/R/news.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


news <-
function(query, package = "R", lib.loc = NULL,
         format = NULL, reader = NULL, db = NULL)
{
    if(is.null(db)) {
        db <- if(package == "R")
            tools:::.build_news_db_from_R_NEWS_Rd()
        else
            tools:::.build_news_db(package, lib.loc, format, reader)
    }
    if(is.null(db))
        return(invisible())

    attr(db, "package") <- package

    ## Is there a way to directly call/use subset.data.frame?
    ## E.g.,
    ##   subset(db, query)
    ## does not work.
    if(missing(query))
        return(db)

    ## For queries we really need to force Version to package_version
    ## and Date to Date ...
    ## This is tricky because we do not necessarily have valid package
    ## versions (e.g., R NEWS has "2.8.1 patched") or could have the
    ## version info missing (and package_version() does not like NAs).

    has_bad_attr <-
        !is.null(bad <- attr(db, "bad")) && (length(bad) == NROW(db))

    ## Manipulate fields for querying (but return the original ones).
    db1 <- db
    ## Canonicalize version entries which *start* with a valid numeric
    ## version, i.e., drop things like " patched".
    version <- db$Version
    pos <- regexpr(sprintf("^%s",
                           .standard_regexps()$valid_numeric_version),
                   version)
    if(any(ind <- (pos > -1L)))
        version[ind] <-
            substring(version[ind], 1L, attr(pos, "match.length")[ind])
    db1$Version <- numeric_version(version, strict = FALSE)
    db1$Date <- as.Date(db$Date)

    r <- eval(substitute(query), db1, parent.frame())
    ## Do something if this is not logical ...
    if(is.null(r))
        return(db)
    else if(!is.logical(r) || length(r) != length(version))
        stop("invalid query")
    r <- r & !is.na(r)
    if(has_bad_attr)
        structure(db[r, ], bad = bad[r])
    else
        db[r, ]
}

format.news_db <-
function(x, ...)
{
    if(inherits(x, "news_db_from_Rd") ||
       (!(is.null(bad <- attr(x, "bad")))
        && (length(bad) == NROW(x))
        && all(!bad))) {

        ## Format news in the preferred input format:
        ##   Changes in $VERSION [($DATE)]:
        ##   [$CATEGORY$]
        ##   indented/formatted bullet list of $TEXT entries.
        ## <FIXME>
        ## Add support for DATE.
        ## </FIXME>

        vchunks <- split(x, x$Version)
        ## Re-order according to decreasing version.
        ## R NEWS has invalid "versions" such as "R-devel" and
        ## "2.4.1 patched".  We can remap the latter (to e.g. 2.4.1.1)
        ## and need to ensure the former come first.
        vstrings <- names(vchunks)
        ind <- vstrings != "R-devel"
        pos <- c(which(!ind),
                 which(ind)[order(as.numeric_version(sub(" *patched", ".1",
                                                         vstrings[ind])),
                                  decreasing = TRUE)])
        vchunks <- vchunks[pos]
	if(length(vchunks)) {
            dates <- sapply(vchunks, function(v) v$Date[1L])
            vstrings <- names(vchunks)
            ind <- vstrings != "R-devel"
            vstrings[ind] <- sprintf("version %s", vstrings[ind])
            vheaders <-
                sprintf("Changes in %s%s:",
                        vstrings,
                        ifelse(is.na(dates), "",
                               sprintf(" (%s)", dates)))
        } else vheaders <- character()

        format_items <- function(x)
            paste0("    o   ", gsub("\n", "\n\t", x$Text))
        format_vchunk <- function(vchunk) {
            if(all(!is.na(category <- vchunk$Category)
                   & nzchar(category))) {
                ## need to preserve order of headings.
                cchunks <-
                    split(vchunk,
                          factor(category, levels = unique(category)))
                Map(c, names(cchunks), sapply(cchunks, format_items),
                    USE.NAMES = FALSE)
            } else {
                format_items(vchunk)
            }
        }

        Map(c, vheaders, lapply(vchunks, format_vchunk),
            USE.NAMES = FALSE)
    } else {
        ## Simple and ugly.
        ## Drop all-NA variables.
        apply(as.matrix(x),
              1L,
              function(e)
              paste(formatDL(e[!is.na(e)], style = "list"),
                    collapse = "\n"))
    }
}

print.news_db <-
function(x, ...)
{
    writeLines(paste(unlist(format(x, ...)), collapse = "\n\n"))
    invisible(x)
}
#  File src/library/utils/R/object.size.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

object.size <- function(x)
    structure(.Call(C_objectSize, x), class = "object_size")

print.object_size <-
    function(x, quote = FALSE, units = "b", ...)
{
    units <- match.arg(units, c("b", "auto", "Kb", "Mb", "Gb",
                                "B", "KB", "MB", "GB"))
    if (units == "auto") {
        if (x >= 1024^3) units <- "Gb"
        else if (x >= 1024^2) units <- "Mb"
        else if (x >= 1024) units <- "Kb"
        else units <- "b"
    }
    y <- switch(units,
                "b" =, "B" = paste(x, "bytes"),
                "Kb" =, "KB" = paste(round(x/1024, 1L), "Kb"),
                "Mb" =, "MB" = paste(round(x/1024^2, 1L), "Mb"),
                "Gb" =, "GB" = paste(round(x/1024^3, 1L), "Gb")
                )
    if(quote) print.default(y, ...) else cat(y, "\n", sep = "")
    invisible(x)
}
#  File src/library/utils/R/objects.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

## findGeneric(fname) :  is 'fname' the name of an S3 generic ?
##			[unexported function used only in this file]
findGeneric <-
function(fname, envir)
{
    if(!exists(fname, mode = "function", envir = envir)) return("")
    f <- get(fname, mode = "function", envir = envir)
    ## FIXME? In the first case, e.g. 'methods(qr)', we are very inefficient:
    ##  inside methods() we transform the 'qr' function object into a character,
    ##  whereas here, we revert this, searching around unnecessarily
    ##
    if(.isMethodsDispatchOn() && methods::is(f, "genericFunction")) {
	## maybe an S3 generic was turned into the S4 default
	## Try to find it, otherwise warn :
	fMethsEnv <- methods::getMethodsForDispatch(f)
	r <- lapply(grep("^ANY\\b", ls(envir = fMethsEnv), value=TRUE),
		    get, envir = fMethsEnv)
	if(any(ddm <- unlist(lapply(r, class)) == "derivedDefaultMethod"))
	    f <- r[ddm][[1]]@.Data
	else
	    warning(gettextf(
	"'%s' is a formal generic function; S3 methods will not likely be found",
			     fname), domain = NA)
    }
    isUMEbrace <- function(e) {
        for (ee in as.list(e[-1L]))
            if (nzchar(res <- isUME(ee))) return(res)
        ""
    }
    isUMEif <- function(e) {
        if (length(e) == 3L) isUME(e[[3L]])
        else {
            if (nzchar(res <- isUME(e[[3L]]))) res
            else if (nzchar(res <- isUME(e[[4L]]))) res
            else ""
        }
    }
    isUME <- function(e) { ## is it an "UseMethod() calling function" ?
        if (is.call(e) && (is.name(e[[1L]]) || is.character(e[[1L]]))) {
            switch(as.character(e[[1L]]),
                   UseMethod = as.character(e[[2L]]),
                   "{" = isUMEbrace(e),
                   "if" = isUMEif(e),
                   "")
        } else ""
    }
    isUME(body(f))
}

getKnownS3generics <-
function()
    c(names(.knownS3Generics), tools:::.get_internal_S3_generics())

methods <-
function(generic.function, class)
{
    rbindSome <- function(df, nms, msg) {
        ## rbind.data.frame() -- dropping rows with duplicated names
        nms <- unique(nms)
        n2 <- length(nms)
        dnew <- data.frame(visible = rep.int(FALSE, n2),
                           from    = rep.int(msg,   n2),
                           row.names = nms)
        n <- nrow(df)
        if(n == 0L) return(dnew)
        ## else
        keep <- !duplicated(c(rownames(df), rownames(dnew)))
        rbind(df  [keep[1L:n] , ],
              dnew[keep[(n+1L):(n+n2)] , ])
    }

    S3MethodsStopList <- tools:::.make_S3_methods_stop_list(NULL)
    knownGenerics <- getKnownS3generics()
    sp <- search()
    an <- lapply(seq_along(sp), ls)
    names(an) <- sp
    an <- unlist(an)
    an <- an[!duplicated(an)] # removed masked objects, *keep* names
    names(an) <- sub("[0-9]*$", "", names(an))
    info <- data.frame(visible = rep.int(TRUE, length(an)),
                       from = names(an),
                       row.names = an)
    if (!missing(generic.function)) {
	if (!is.character(generic.function))
	    generic.function <- deparse(substitute(generic.function))
        ## else
        if(!exists(generic.function, mode = "function",
                   envir = parent.frame()) &&
           !any(generic.function == c("Math", "Ops", "Complex", "Summary")))
            stop(gettextf("no function '%s' is visible", generic.function),
                 domain = NA)
        if(!any(generic.function == knownGenerics)) {
            truegf <- findGeneric(generic.function, parent.frame())
            if(truegf == "")
                warning(gettextf("function '%s' appears not to be generic",
                                 generic.function), domain = NA)
            else if(truegf != generic.function) {
                warning(gettextf("generic function '%s' dispatches methods for generic '%s'",
                        generic.function, truegf), domain = NA)
                generic.function <- truegf
            }
        }
	name <- paste0("^", generic.function, ".")
        name <- gsub("([.[$+*])", "\\\\\\1",name)
        info <- info[grep(name, row.names(info)), ]
        info <- info[! row.names(info) %in% S3MethodsStopList, ]
        ## check that these are all functions
        ## might be none at this point
        if(nrow(info)) {
            keep <- sapply(row.names(info),
                           function(nm) exists(nm, mode="function"))
            info <- info[keep, ]
        }

        ## also look for registered methods from namespaces
        ## we assume that only functions get registered.
        defenv <- if(!is.na(w <- .knownS3Generics[generic.function]))
            asNamespace(w)
        else {
            genfun <- get(generic.function, mode = "function",
                          envir = parent.frame())
            if(.isMethodsDispatchOn() && methods::is(genfun, "genericFunction"))
                genfun <- methods::finalDefaultMethod(genfun@default)
            if (typeof(genfun) == "closure") environment(genfun)
            else .BaseNamespaceEnv
        }
        S3reg <- ls(get(".__S3MethodsTable__.", envir = defenv),
                    pattern = name)
        if(length(S3reg))
            info <- rbindSome(info, S3reg, msg =
                              paste("registered S3method for",
                                    generic.function))
        ## both all() and all.equal() are generic, so
        if(generic.function == "all")
            info <- info[-grep("^all\\.equal", row.names(info)), ]
    }
    else if (!missing(class)) {
	if (!is.character(class))
	    class <- paste(deparse(substitute(class)))
	name <- paste0(".", class, "$")
        name <- gsub("([.[])", "\\\\\\1", name)
        info <- info[grep(name, row.names(info)), ]
        info <- info[! row.names(info) %in% S3MethodsStopList, ]

        if(nrow(info)) {
            ## check if we can find a generic matching the name
            possible.generics <- gsub(name, "", row.names(info))
            keep <- sapply(possible.generics, function(nm) {
                if(nm %in% knownGenerics) return(TRUE)
                where <- find(nm, mode = "function")
                if(!length(where)) return(FALSE)
                any(sapply(where, function(w)
                           nzchar(findGeneric(nm, envir=as.environment(w)))))
            })
            info <- info[keep, ]
        }

        ## also look for registered methods in loaded namespaces.
        ## These should only be registered in environments containing
        ## the corresponding generic, so we don't check again.
        ## Note that the generic will not necessarily be visible,
        ## as the package may not be loaded.
        S3reg <- unlist(lapply(loadedNamespaces(), function(i) ls(get(".__S3MethodsTable__.", envir = asNamespace(i)), pattern = name)))
        ## now methods like print.summary.aov will be picked up,
        ## so we do look for such mismatches.
        if(length(S3reg))
            S3reg <- S3reg[sapply(gsub(name, "", S3reg), exists)]
        if(length(S3reg))
            info <- rbindSome(info, S3reg, msg = "registered S3method")
    }
    else stop("must supply 'generic.function' or 'class'")

    info <- info[sort.list(row.names(info)), ]
    res <- row.names(info)
    class(res) <- "MethodsFunction"
    attr(res, "info") <- info
    res
}

print.MethodsFunction <-
function(x, ...)
{
    visible <- attr(x, "info")[["visible"]]
    if(length(x)) {
	print(paste0(x, ifelse(visible, "", "*")), quote=FALSE, ...)
        if(any(!visible))
            cat("\n", "   ",
                "Non-visible functions are asterisked", "\n", sep = "")
    } else cat("no methods were found\n")
    invisible(x)
}


getS3method <-
function(f, class, optional = FALSE)
{
    if(!any(f == getKnownS3generics())) {
        truegf <- findGeneric(f, parent.frame())
        if(nzchar(truegf)) f <- truegf
        else {
            if(optional) return(NULL)
            else stop(gettextf("no function '%s' could be found", f), domain = NA)
        }
    }
    method <- paste(f, class, sep=".")
    if(exists(method, mode = "function", envir = parent.frame()))
        return(get(method, mode = "function", envir = parent.frame()))
    ## also look for registered method in namespaces
    defenv <- if(!is.na(w <- .knownS3Generics[f])) asNamespace(w)
    else if(f %in% tools:::.get_internal_S3_generics()) .BaseNamespaceEnv
    else {
        genfun <- get(f, mode="function", envir = parent.frame())
        if(.isMethodsDispatchOn() && methods::is(genfun, "genericFunction"))
            ## assumes the default method is the S3 generic function
            genfun <- methods::selectMethod(genfun, "ANY")
        if (typeof(genfun) == "closure") environment(genfun)
        else .BaseNamespaceEnv
    }
    S3Table <- get(".__S3MethodsTable__.", envir = defenv)
    if(exists(method, envir = S3Table, inherits = FALSE))
        return(get(method, envir = S3Table))
    if(optional) NULL else stop(gettextf("S3 method '%s' not found", method),
                                domain = NA)
}

getFromNamespace <-
function(x, ns, pos = -1, envir = as.environment(pos))
{
    if(missing(ns)) {
        nm <- attr(envir, "name", exact = TRUE)
        if(is.null(nm) || substring(nm, 1L, 8L) != "package:")
            stop("environment specified is not a package")
        ns <- asNamespace(substring(nm, 9L))
    } else ns <- asNamespace(ns)
    get(x, envir = ns, inherits = FALSE)
}

assignInMyNamespace <-
function(x, value)
{
    f <- sys.function(-1)
    ns <- environment(f)
    ## deal with subclasses of "function"
    ## that may insert an environment in front of the namespace
    if(isS4(f))
        while(!isNamespace(ns))
            ns <- parent.env(ns)
    if(bindingIsLocked(x, ns)) {
        unlockBinding(x, ns)
        assign(x, value, envir = ns, inherits = FALSE)
        w <- options("warn")
        on.exit(options(w))
        options(warn = -1)
        lockBinding(x, ns)
    } else assign(x, value, envir = ns, inherits = FALSE)
    if(!isBaseNamespace(ns)) {
        ## now look for possible copy as a registered S3 method
        S3 <- getNamespaceInfo(ns, "S3methods")
        if(!length(S3)) return(invisible(NULL))
        S3names <- S3[, 3L]
        if(x %in% S3names) {
            i <- match(x, S3names)
            genfun <- get(S3[i, 1L], mode = "function", envir = parent.frame())
            if(.isMethodsDispatchOn() && methods::is(genfun, "genericFunction"))
                genfun <- methods::slot(genfun, "default")@methods$ANY
            defenv <- if (typeof(genfun) == "closure") environment(genfun)
            else .BaseNamespaceEnv
            S3Table <- get(".__S3MethodsTable__.", envir = defenv)
            remappedName <- paste(S3[i, 1L], S3[i, 2L], sep = ".")
            if(exists(remappedName, envir = S3Table, inherits = FALSE))
                assign(remappedName, value, S3Table)
        }
    }
    invisible(NULL)
}

assignInNamespace <-
function(x, value, ns, pos = -1, envir = as.environment(pos))
{
    nf <- sys.nframe()
    if(missing(ns)) {
        nm <- attr(envir, "name", exact = TRUE)
        if(is.null(nm) || substring(nm, 1L, 8L) != "package:")
            stop("environment specified is not a package")
        ns <- asNamespace(substring(nm, 9L))
    } else ns <- asNamespace(ns)
    if (nf > 1L) {
        if(getNamespaceName(ns) %in% tools:::.get_standard_package_names()$base)
            stop("locked binding of ", sQuote(x), " cannot be changed",
                 domain = NA)
    }
    if(bindingIsLocked(x, ns)) {
        in_load <- Sys.getenv("_R_NS_LOAD_")
        if (nzchar(in_load)) {
            ns_name <- getNamespaceName(ns)
            if(in_load != ns_name) {
                msg <-
                    gettextf("changing locked binding for %s in %s whilst loading %s",
                             sQuote(x), sQuote(ns_name), sQuote(in_load))
                if (! in_load %in% c("Matrix", "SparseM"))
                    warning(msg, call. = FALSE, domain = NA, immediate. = TRUE)
            }
        } else if (nzchar(Sys.getenv("_R_WARN_ON_LOCKED_BINDINGS_"))) {
            ns_name <- getNamespaceName(ns)
            warning(gettextf("changing locked binding for %s in %s",
                             sQuote(x), sQuote(ns_name)),
                    call. = FALSE, domain = NA, immediate. = TRUE)
        }
        unlockBinding(x, ns)
        assign(x, value, envir = ns, inherits = FALSE)
        w <- options("warn")
        on.exit(options(w))
        options(warn = -1)
        lockBinding(x, ns)
    } else {
        assign(x, value, envir = ns, inherits = FALSE)
    }
    if(!isBaseNamespace(ns)) {
        ## now look for possible copy as a registered S3 method
        S3 <- getNamespaceInfo(ns, "S3methods")
        if(!length(S3)) return(invisible(NULL))
        S3names <- S3[, 3L]
        if(x %in% S3names) {
            i <- match(x, S3names)
            genfun <- get(S3[i, 1L], mode = "function", envir = parent.frame())
            if(.isMethodsDispatchOn() && methods::is(genfun, "genericFunction"))
                genfun <- methods::slot(genfun, "default")@methods$ANY
            defenv <- if (typeof(genfun) == "closure") environment(genfun)
            else .BaseNamespaceEnv
            S3Table <- get(".__S3MethodsTable__.", envir = defenv)
            remappedName <- paste(S3[i, 1L], S3[i, 2L], sep = ".")
            if(exists(remappedName, envir = S3Table, inherits = FALSE))
                assign(remappedName, value, S3Table)
        }
    }
    invisible(NULL)
}

fixInNamespace <-
function(x, ns, pos = -1, envir = as.environment(pos), ...)
{
    subx <- substitute(x)
    if (is.name(subx))
        subx <- deparse(subx)
    if (!is.character(subx) || length(subx) != 1L)
        stop("'fixInNamespace' requires a name")
    if(missing(ns)) {
        nm <- attr(envir, "name", exact = TRUE)
        if(is.null(nm) || substring(nm, 1L, 8L) != "package:")
            stop("environment specified is not a package")
        ns <- asNamespace(substring(nm, 9L))
    } else ns <- asNamespace(ns)
    x <- edit(get(subx, envir = ns, inherits = FALSE), ...)
    assignInNamespace(subx, x, ns)
}

getAnywhere <-
function(x)
{
    if(tryCatch(!is.character(x), error = function(e) TRUE))
        x <- as.character(substitute(x))
    objs <- list(); where <- character(); visible <- logical()
    ## first look on search path
    if(length(pos <- find(x, numeric = TRUE))) {
        objs <- lapply(pos, function(pos, x) get(x, pos=pos), x=x)
        where <- names(pos)
        visible <- rep.int(TRUE, length(pos))
    }
    ## next look for methods: a.b.c.d could be a method for a or a.b or a.b.c
    if(length(grep(".", x, fixed=TRUE))) {
        np <- length(parts <- strsplit(x, ".", fixed=TRUE)[[1L]])
        for(i in 2:np) {
            gen <- paste(parts[1L:(i-1)], collapse=".")
            cl <- paste(parts[i:np], collapse=".")
            if (gen == "" || cl == "") next
            ## want to evaluate this in the parent, or the utils namespace
            ## gets priority.
            Call <- substitute(getS3method(gen, cl, TRUE), list(gen = gen, cl = cl))
            f <- eval.parent(Call)
            ## Now try to fathom out where it is from.
            ## f might be a special, not a closure, and not have an environment,
            if(!is.null(f) && !is.null(environment(f))) {
                ev <- topenv(environment(f), baseenv())
                nmev <- if(isNamespace(ev)) getNamespaceName(ev) else NULL
                objs <- c(objs, f)
                msg <- paste("registered S3 method for", gen)
                if(!is.null(nmev))
                    msg <- paste(msg, "from namespace", nmev)
                where <- c(where, msg)
                visible <- c(visible, FALSE)
            }
        }
    }
    ## now look in loaded namespaces
    for(i in loadedNamespaces()) {
        ns <- asNamespace(i)
        if(exists(x, envir = ns, inherits = FALSE)) {
            f <- get(x, envir = ns, inherits = FALSE)
            objs <- c(objs, f)
            where <- c(where, paste("namespace", i, sep=":"))
            visible <- c(visible, FALSE)
        }
    }
    # now check for duplicates
    ln <- length(objs)
    dups <- rep.int(FALSE, ln)
    if(ln > 1L)
        for(i in 2L:ln)
            for(j in 1L:(i-1L))
                if(identical(objs[[i]], objs[[j]],
                             ignore.environment = TRUE)) {
                    dups[i] <- TRUE
                    break
                }
    res <- list(name=x, objs=objs, where=where, visible=visible, dups=dups)
    class(res) <- "getAnywhere"
    res
}

print.getAnywhere <-
function(x, ...)
{
    n <- sum(!x$dups)
    if(n == 0L) {
        cat("no object named", sQuote(x$name), "was found\n")
    } else if (n == 1L) {
        cat("A single object matching", sQuote(x$name), "was found\n")
        cat("It was found in the following places\n")
	cat(paste0("  ", x$where), sep="\n")
        cat("with value\n\n")
        print(x$objs[[1L]])
    } else {
        cat(n, "differing objects matching", sQuote(x$name),
            "were found\n")
        cat("in the following places\n")
        cat(paste0("  ", x$where), sep="\n")
        cat("Use [] to view one of them\n")
    }
    invisible(x)
}

`[.getAnywhere` <-
function(x, i)
{
    if(!is.numeric(i)) stop("only numeric indices can be used")
    if(length(i) == 1L) x$objs[[i]]
    else x$objs[i]
}

argsAnywhere <-
function(x)
{
    if(tryCatch(!is.character(x), error = function(e) TRUE))
        x <- as.character(substitute(x))
    fs <- getAnywhere(x)
    if (sum(!fs$dups) == 0L)
        return(NULL)
    if (sum(!fs$dups) > 1L)
        sapply(fs$objs[!fs$dups],
               function(f) if (is.function(f)) args(f))
    else args(fs$objs[[1L]])
}
#  File src/library/utils/R/package.skeleton.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

package.skeleton <-
    function(name = "anRpackage", list = character(), environment = .GlobalEnv,
	     path = ".", force = FALSE, namespace = TRUE,
             code_files = character())
{
    safe.dir.create <- function(path)
    {
	dirTest <- function(x) !is.na(isdir <- file.info(x)$isdir) & isdir
	if(!dirTest(path) && !dir.create(path))
	    stop(gettextf("cannot create directory '%s'", path), domain = NA)
    }

    if(!is.character(code_files))
        stop("'code_files' must be a character vector")
    use_code_files <- length(code_files) > 0L

    envIsMissing <- missing(environment) # before R clobbers this information

    if(missing(list)) {
        if(use_code_files) {
            environment <- new.env(hash = TRUE)
            methods::setPackageName(name, environment)
            for(cf in code_files)
                sys.source(cf, envir = environment)
        }
	## all.names: crucial for metadata
	list <- ls(environment, all.names=TRUE)
    }
    if(!is.character(list))
	stop("'list' must be a character vector naming R objects")
    if(use_code_files || !envIsMissing) {
        classesList <- getClasses(environment)
        classes0 <- .fixPackageFileNames(classesList)
        names(classes0) <- classesList
        methodsList <- getGenerics(environment)
        methods0 <- .fixPackageFileNames(methodsList)
        names(methods0) <- methodsList
    }
    else { # nobody should  specify classes or methods as object names!
        classesList <- methodsList <- character()
    }
    usingS4 <- length(classesList) > 0L || length(methodsList) > 0L

    ## we need to test in the C locale
    curLocale <- Sys.getlocale("LC_CTYPE")
    on.exit(Sys.setlocale("LC_CTYPE", curLocale), add = TRUE)
    if(Sys.setlocale("LC_CTYPE", "C") != "C")
        warning("cannot turn off locale-specific chars via LC_CTYPE",
                domain = NA)

    have <- unlist(lapply(list, exists, envir = environment))
    if(any(!have))
        warning(sprintf(ngettext(sum(!have),
                                 "object '%s' not found",
                                 "objects '%s' not found"),
                        paste(sQuote(list[!have]), collapse=", ")),
                domain = NA)
    list <- list[have]
    if(!length(list))
	stop("no R objects specified or available")

    message("Creating directories ...", domain = NA)
    ## Make the directories
    dir <- file.path(path, name)
    if(file.exists(dir) && !force)
	stop(gettextf("directory '%s' already exists", dir), domain = NA)

    safe.dir.create(dir)
    safe.dir.create(code_dir <- file.path(dir, "R"))
    safe.dir.create(docs_dir <- file.path(dir, "man"))
    safe.dir.create(data_dir <- file.path(dir, "data"))

    ## DESCRIPTION
    message("Creating DESCRIPTION ...", domain = NA)
    description <- file(file.path(dir, "DESCRIPTION"), "wt")
    cat("Package: ", name, "\n",
	"Type: Package\n",
	"Title: What the package does (short line)\n",
	"Version: 1.0\n",
	"Date: ", format(Sys.time(), format="%Y-%m-%d"), "\n",
	"Author: Who wrote it\n",
	"Maintainer: Who to complain to <yourfault@somewhere.net>\n",
	"Description: More about what it does (maybe more than one line)\n",
	"License: What license is it under?\n",
	if(usingS4) "Depends: methods\n",
	file = description, sep = "")
    close(description)

    if(!missing(namespace))
	warning("From R 2.14.0 on, every package gets a NAMESPACE.",
		" Argument 'namespace' is deprecated.", domain = NA)
    ## NAMESPACE
    ## <NOTE>
    ## For the time being, we export all non-internal objects using the pattern
    ## of names beginning with alpha.  All S4 methods and classes are exported.
    ## S3 methods will be exported if the function's name would be exported.
    ## </NOTE>
    message("Creating NAMESPACE ...", domain = NA)
    out <- file(file.path(dir, "NAMESPACE"), "wt")
    writeLines("exportPattern(\"^[[:alpha:]]+\")", out)
    if(length(methodsList)) {
	cat("exportMethods(\n    ", file = out)
	cat(paste0('"', methodsList, '"', collapse = ",\n    "), "\n)\n", file = out)
    }
    if(length(classesList)) {
	cat("exportClasses(\n    ", file = out)
	cat(paste0('"', classesList, '"', collapse = ",\n     "), "\n)\n", file = out)
    }
    close(out)

    ## Read-and-delete-me
    message("Creating Read-and-delete-me ...", domain = NA)
    out <- file(file.path(dir, "Read-and-delete-me"), "wt")
    msg <-
        c("* Edit the help file skeletons in 'man', possibly combining help files for multiple functions.",
          "* Edit the exports in 'NAMESPACE', and add necessary imports.",
          "* Put any C/C++/Fortran code in 'src'.",
          "* If you have compiled code, add a useDynLib() directive to 'NAMESPACE'.",
          "* Run R CMD build to build the package tarball.",
          "* Run R CMD check to check the package tarball.",
          "",
          "Read \"Writing R Extensions\" for more information.")
    writeLines(strwrap(msg, exdent = 2), out)
    close(out)

    internalObjInds <- grep("^\\.", list)
    internalObjs <- list[internalObjInds]
    if(length(internalObjInds))
	list <- list[-internalObjInds]

    list0 <- .fixPackageFileNames(list)
    names(list0) <- list

    ## Dump the items in 'data' or 'R'
    if(!use_code_files) {
        message("Saving functions and data ...", domain = NA)
        if(length(internalObjInds))
            dump(internalObjs,
                 file = file.path(code_dir, sprintf("%s-internal.R", name)),
                 envir = environment)
        for(item in list){
            objItem <- get(item, envir = environment)
            if(is.function(objItem))  {
                if(isS4(objItem))
                    stop(gettextf("generic functions and other S4 objects (e.g., '%s') cannot be dumped; use the 'code_files' argument", item), domain = NA)
                dump(item,
                     file = file.path(code_dir, sprintf("%s.R", list0[item])),
                     envir = environment)
            }
            else       # we cannot guarantee this is a valid file name
                try(save(list = item, envir = environment,
                         file = file.path(data_dir, sprintf("%s.rda", item))))
        }
    } else {
        message("Copying code files ...", domain = NA)
        file.copy(code_files, code_dir)
        ## Only "abc.R"-like files are really ok:
	R_files <- tools::list_files_with_type(code_dir, "code",
					       full.names = FALSE,
					       OS_subdirs = "")
        code_files <- basename(code_files)
	wrong <- code_files[is.na(match(code_files, R_files))]
	if(length(wrong)) {
	    warning("Invalid file name(s) for R code in ", code_dir,":\n",
		    strwrap(paste(sQuote(wrong), collapse = ", "), indent=2),
		    "\n are now renamed to 'z<name>.R'", domain = NA)
	    file.rename(from = file.path(code_dir, wrong),
			to = file.path(code_dir,
			paste0("z", sub("(\\.[^.]*)?$", ".R", wrong))))
        }
    }

    ## Make help file skeletons in 'man'
    message("Making help files ...", domain = NA)
    ## Suppress partially inappropriate messages from prompt().
    yy <- try(suppressMessages({
	promptPackage(name,
		      filename =
		      file.path(docs_dir,
				sprintf("%s-package.Rd", name)),
		      lib.loc = path)
	sapply(list,
	       function(item) {
		   prompt(get(item, envir = environment),
			  name = item,
			  filename =
			  file.path(docs_dir,
				    sprintf("%s.Rd", list0[item])))
	       })
	sapply(classesList,
	       function(item) {
		   methods::promptClass(item,
					filename =
					file.path(docs_dir,
						  sprintf("%s-class.Rd", classes0[item])),
					where = environment)
	       })
	sapply(methodsList,
	       function(item) {
		   methods::promptMethods(item,
					  filename =
					  file.path(docs_dir,
						    sprintf("%s-methods.Rd", methods0[item])),
					  findMethods(item, where = environment))
	       })
    }))
    ## don't document generic functions from other packages
    for(item in methodsList) {
        if(exists(item, envir = environment, inherits = FALSE)) {
            ff <- get(item, envir = environment)
            if(is(ff, "genericFunction") && !identical(ff@package, name)) # don't document
                file.remove(file.path(docs_dir, sprintf("%s.Rd", list0[item])))
        }
    }
    if(inherits(yy, "try-error"))
	stop(yy)

    ## Now we may have created an empty data or R directory
    if(length(list.files(code_dir)) == 0L)
        unlink(code_dir, recursive = TRUE)
    if(length(list.files(data_dir)) == 0L)
        unlink(data_dir, recursive = TRUE)

    message("Done.", domain = NA)
    message(sprintf("Further steps are described in '%s'.",
                     file.path(dir, "Read-and-delete-me")),
            domain = NA)
}

.fixPackageFileNames <- function(list) {
        ## Some object names may not be valid file names, especially
        ## replacement function names.  And if we start changing them
        ## they may collide.
        ## <NOTE>
        ## If we use given code files, we could still check whether
        ## these file are valid across platforms ...
        ## </NOTE>
        list <- as.character(list) # remove S4 class if any, to add names() later
        if(length(list) == 0L) return(list)
        list0 <- gsub("[[:cntrl:]\"*/:<>?\\|]", "_", list)
        wrong <- grep("^(con|prn|aux|clock\\$|nul|lpt[1-3]|com[1-4])(\\..*|)$",
                      list0)
        if(length(wrong))
            list0[wrong] <- paste0("zz", list0[wrong])
        ## using grep was wrong, as could give -integer(0)
        ok <- grepl("^[[:alnum:]]", list0)
        if(any(!ok))
            list0[!ok] <- paste0("z", list0[!ok])
        ## now on Mac/Windows lower/uppercase will collide too
        list1 <- tolower(list0)
        list2 <- make.unique(list1, sep = "_")
        changed <- (list2 != list1)
        list0[changed] <- list2[changed]
        list0
}
#  File src/library/utils/R/packageStatus.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

packageStatus <- function(lib.loc = NULL, repositories = NULL, method,
                          type = getOption("pkgType"))
{
    newestVersion <- function(x)
    {
        vers <- package_version(x)
	max <- vers[1L]
        for (i in seq_along(vers)) if (max < vers[i]) max <- vers[i]
	which.max(vers == max)
    }

    if(is.null(lib.loc))
        lib.loc <- .libPaths()
    if(is.null(repositories))
        repositories <- contrib.url(getOption("repos"), type = type)

    ## convert character matrices to dataframes
    char2df <- function(x)
    {
        y <- list()
        for(k in 1L:ncol(x)) y[[k]] <- x[,k]
        attr(y, "names") <- colnames(x)
        attr(y, "row.names") <- make.unique(y[[1L]])
        class(y) <- "data.frame"
        y
    }

    y <- char2df(installed.packages(lib.loc = lib.loc))
    y[, "Status"] <- "ok"

    z <- available.packages(repositories, method)
    ## only consider the newest version of each package
    ## in the first repository where it appears
    ztab <- table(z[,"Package"])
    for(pkg in names(ztab)[ztab>1]){
        zrow <- which(z[,"Package"]==pkg)
        znewest <- newestVersion(z[zrow,"Version"])
        ## and now exclude everything but the newest
        z <- z[-zrow[-znewest],]
    }

    z <- cbind(z, Status = "not installed")
    z[z[,"Package"] %in% y$Package, "Status"] <- "installed"

    z <- char2df(z)
    attr(z, "row.names") <- z$Package

    for(k in 1L:nrow(y)){
        pkg <- y[k, "Package"]
        if(pkg %in% z$Package) {
            if(package_version(y[k, "Version"]) <
               package_version(z[pkg, "Version"])) {
                y[k, "Status"] <- "upgrade"
            }
        } else {
            if(!(y[k, "Priority"] %in% "base")) y[k, "Status"] <- "unavailable"
        }
    }

    y$LibPath <- factor(y$LibPath, levels=lib.loc)
    y$Status <- factor(y$Status, levels=c("ok", "upgrade", "unavailable"))
    z$Repository <- factor(z$Repository, levels=repositories)
    z$Status <- factor(z$Status, levels=c("installed", "not installed"))

    retval <- list(inst=y, avail=z)
    class(retval) <- "packageStatus"
    retval
}

summary.packageStatus <- function(object, ...)
{
    Libs <- levels(object$inst$LibPath)
    Repos <- levels(object$avail$Repository)

    byLib <- split(object$inst, object$inst$LibPath)
    Libs <- lapply(split(object$inst, object$inst$LibPath),
                   function(x) tapply(x$Package, x$Status,
                                      function(x) sort(as.character(x))))
    Repos <- lapply(split(object$avail, object$avail$Repository),
                    function(x) tapply(x$Package, x$Status,
                                       function(x) sort(as.character(x))))
    object$Libs <- Libs
    object$Repos <- Repos
    class(object) <- c("summary.packageStatus", "packageStatus")
    object
}

print.summary.packageStatus <- function(x, ...)
{
    cat("\nInstalled packages:\n")
    cat(  "-------------------\n")
    for(k in seq_along(x$Libs)) {
        cat("\n*** Library ", names(x$Libs)[k], "\n", sep = "")
	print(x$Libs[[k]], ...)
    }
    cat("\n\nAvailable packages:\n")
    cat(    "-------------------\n")
    cat("(each package appears only once)\n")
    for(k in seq_along(x$Repos)){
        cat("\n*** Repository ", names(x$Repos)[k], "\n", sep = "")
	print(x$Repos[[k]], ...)
    }
    invisible(x)
}

print.packageStatus <- function(x, ...)
{
    cat("Number of installed packages:\n")
    print(table(x$inst$LibPath, x$inst$Status), ...)

    cat("\nNumber of available packages (each package counted only once):\n")
    print(table(x$avail$Repository, x$avail$Status), ...)
    invisible(x)
}

update.packageStatus <-
    function(object, lib.loc=levels(object$inst$LibPath),
             repositories=levels(object$avail$Repository),
             ...)
{
    packageStatus(lib.loc=lib.loc, repositories=repositories)
}


upgrade <- function(object, ...)
    UseMethod("upgrade")

upgrade.packageStatus <- function(object, ask=TRUE, ...)
{
    update <- NULL
    old <- which(object$inst$Status == "upgrade")
    if(length(old) == 0L) {
        cat("Nothing to do!\n")
        return(invisible())
    }

    askprint <- function(x)
        write.table(x, row.names = FALSE, col.names = FALSE, quote = FALSE,
                    sep = " at ")

    haveasked <- character()
    if(ask) {
        for(k in old) {
            pkg <-  object$inst[k, "Package"]
            tmpstring <- paste(pkg, as.character(object$inst[k, "LibPath"]))
            if(tmpstring %in% haveasked) next
            haveasked <- c(haveasked, tmpstring)
            cat("\n")
            cat(pkg, ":\n")
            askprint(object$inst[k,c("Version", "LibPath")])
            askprint(object$avail[pkg, c("Version", "Repository")])
            answer <- substr(readline("Update (y/N/x)?  "), 1L, 1L)
            if(answer == "c" | answer == "c") {
                cat("cancelled by user\n")
                return(invisible())
            }
            if(answer == "y" | answer == "Y")
                update <-
                    rbind(update,
                          c(pkg, as.character(object$inst[k, "LibPath"]),
                            as.character(object$avail[pkg, "Repository"])))
        }
    } else {
        pkgs <- object$inst[ ,"Package"]
        update <- cbind(pkgs, as.character(object$inst[ , "LibPath"]),
                        as.character(object$avail[pkgs, "Repository"]))
        update <- update[old, , drop=FALSE]
    }

    if(length(update)) {
        for(repo in unique(update[,3])) {
            ok <- update[, 3] == repo
            install.packages(update[ok, 1], update[ok, 2], contriburl = repo)
        }
    }
}
#  File src/library/utils/R/packages.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

available.packages <-
function(contriburl = contrib.url(getOption("repos"), type), method,
         fields = NULL, type = getOption("pkgType"), filters = NULL)
{
    requiredFields <-
        c(tools:::.get_standard_repository_db_fields(), "File")
    if (is.null(fields))
	fields <- requiredFields
    else {
	stopifnot(is.character(fields))
	fields <- unique(c(requiredFields, fields))
    }

    res <- matrix(NA_character_, 0L, length(fields) + 1L,
		  dimnames = list(NULL, c(fields, "Repository")))

    for(repos in contriburl) {
        localcran <- length(grep("^file:", repos)) > 0L
        if(localcran) {
            ## see note in download.packages
            if(substring(repos, 1L, 8L) == "file:///") {
                tmpf <- paste(substring(repos, 8L), "PACKAGES", sep = "/")
                if(.Platform$OS.type == "windows") {
                    if(length(grep("^/[A-Za-z]:", tmpf)))
                        tmpf <- substring(tmpf, 2L)
                }
            } else {
                tmpf <- paste(substring(repos, 6L), "PACKAGES", sep = "/")
            }
            res0 <- read.dcf(file = tmpf)
            if(length(res0)) rownames(res0) <- res0[, "Package"]
        } else {
            dest <- file.path(tempdir(),
                              paste0("repos_", URLencode(repos, TRUE), ".rds"))
            if(file.exists(dest)) {
                res0 <- readRDS(dest)
            } else {
                tmpf <- tempfile()
                on.exit(unlink(tmpf))
                op <- options("warn")
                options(warn = -1)

                ## Two kinds of errors can happen:  PACKAGES.gz may not exist,
                ## or a junk error page that is not a valid dcf file may be
                ## returned.  Handle both...

                ## This is a binary file
                z <- tryCatch(download.file(url = paste(repos, "PACKAGES.gz", sep = "/"),
                                            destfile = tmpf, method = method,
                                            cacheOK = FALSE, quiet = TRUE, mode = "wb"),
                              error = identity)
		if(!inherits(z, "error"))
		    z <- res0 <- tryCatch(read.dcf(file = tmpf), error = identity)
                if(inherits(z, "error")) {
                    ## read.dcf is going to interpret CRLF as LF, so use
                    ## binary mode to avoid CRCRLF.
                    z <- tryCatch(download.file(url = paste(repos, "PACKAGES", sep = "/"),
                                                destfile = tmpf, method = method,
                                                cacheOK = FALSE, quiet = TRUE,
                                                mode = "wb"),
                                  error = identity)
		    options(op)
		    if(inherits(z, "error")) {
			warning(gettextf("unable to access index for repository %s", repos),
				call. = FALSE, immediate. = TRUE, domain = NA)
			next
		    }
		    res0 <- read.dcf(file = tmpf)
		} else
		    options(op)
                ## Do we want to cache an empty result?
                if(length(res0)) rownames(res0) <- res0[, "Package"]
                saveRDS(res0, dest, compress = TRUE)
                unlink(tmpf)
                on.exit()
            } # end of download vs cached
        } # end of localcran vs online
        if (length(res0)) {
            missingFields <- fields[!(fields %in% colnames(res0))]
            if (length(missingFields)) {
                toadd <- matrix(NA_character_, nrow=nrow(res0),
                                ncol=length(missingFields),
                                dimnames=list(NULL, missingFields))
                res0 <- cbind(res0, toadd)
            }
            if ("Path" %in% colnames(res0)) {
                rp <- rep.int(repos, nrow(res0))
                path <- res0[, "Path"]
                rp[!is.na(path)] <- paste(repos, path[!is.na(path)], sep = "/")
            } else rp <- repos
            res0 <- cbind(res0[, fields, drop = FALSE], Repository = rp)
            res <- rbind(res, res0)
        }
    }

    if(!length(res)) return(res)

    if(is.null(filters)) {
        filters <- getOption("available_packages_filters")
        if(is.null(filters))
            filters <- available_packages_filters_default
    }
    if(is.list(filters)) {
        ## If filters is a list with an add = TRUE element, add the
        ## given filters to the default ones.
        if(identical(filters$add, TRUE)) {
            filters$add <- NULL
            filters <- c(available_packages_filters_default, filters)
        }
    }
    for(f in filters) {
        if(!length(res)) break
        if(is.character(f)) {
            ## Look up the filters db.
            ## Could be nice and allow abbrevs or ignore case.
            f <- available_packages_filters_db[[f[1L]]]
        }
        if(!is.function(f))
            stop("invalid 'filters' argument.")
        res <- f(res)
    }

    res
}

available_packages_filters_default <-
    c("R_version", "OS_type", "subarch", "duplicates")

available_packages_filters_db <- new.env(hash = FALSE) # small

available_packages_filters_db$R_version <-
function(db)
{
    ## Ignore packages which don't fit our version of R.
    depends <- db[, "Depends"]
    depends[is.na(depends)] <- ""
    ## Collect the (versioned) R depends entries.
    x <- lapply(strsplit(sub("^[[:space:]]*", "", depends),
                             "[[:space:]]*,[[:space:]]*"),
                function(s) s[grepl("^R[[:space:]]*\\(", s)])
    lens <- sapply(x, length)
    pos <- which(lens > 0L)
    if(!length(pos)) return(db)
    lens <- lens[pos]
    ## Unlist.
    x <- unlist(x)
    pat <- "^R[[:space:]]*\\(([[<>=!]+)[[:space:]]+(.*)\\)[[:space:]]*"
    ## Extract ops.
    ops <- sub(pat, "\\1", x)
    ## Split target versions accordings to ops.
    v_t <- split(sub(pat, "\\2", x), ops)
    ## Current R version.
    v_c <- getRversion()
    ## Compare current to target grouped by op.
    res <- logical(length(x))
    for(op in names(v_t))
        res[ops == op] <- do.call(op, list(v_c, v_t[[op]]))
    ## And assemble test results according to the rows of db.
    ind <- rep.int(TRUE, NROW(db))
    ind[pos] <- sapply(split(res, rep.int(seq_along(lens), lens)), all)
    db[ind, , drop = FALSE]
}

available_packages_filters_db$OS_type <-
function(db)
{
    ## Ignore packages that do not fit our OS.
    OS_type <- db[, "OS_type"]
    db[is.na(OS_type) | (OS_type == .Platform$OS.type), , drop = FALSE]
}

available_packages_filters_db$subarch <-
function(db)
{
    ## Ignore packages that do not fit our sub-architecture.
    ## Applies only to Mac and Windows binary repositories.
    current <- .Platform$r_arch
    if(!nzchar(current)) return(db)
    archs <- db[, "Archs"]
    if(all(is.na(archs))) return(db)
    OK <- unlist(lapply(archs, function(x) {
        if(is.na(x)) return(TRUE)
        this <- strsplit(x, "[[:space:]]*,[[:space:]]*")[[1L]]
        current %in% this
    }))
    db[OK, , drop = FALSE]
}

available_packages_filters_db$duplicates <-
function(db)
    tools:::.remove_stale_dups(db)

filter_packages_by_depends_predicates <-
function(db, predicate, recursive = TRUE)
{
    ## Could also add a 'which' argument to specify which dependencies
    ## are taken.

    ## Drop all packages for which any (recursive) dependency does not
    ## satisfy the given predicate (implemented as a function computing
    ## TRUE or FALSE for each rows of the package db).

    ## Somewhat tricky because there may be depends missing from the db,
    ## which are taken not to satisfy the predicate unless they are
    ## standard packages.

    ## Determine all depends missing from the db.
    db1 <- data.frame(Package = db[, "Package"],
                      stringsAsFactors = FALSE)
    fields <- c("Depends", "Imports", "LinkingTo")
    for(f in fields)
        db1[[f]] <-
            lapply(db[, f], tools:::.extract_dependency_package_names)
    all_packages <- unique(unlist(db1[fields], use.names = FALSE))
    bad_packages <-
        all_packages[is.na(match(all_packages, db1$Package))]
    ## Drop the standard packages from these.
    bad_packages <-
        setdiff(bad_packages,
                unlist(tools:::.get_standard_package_names()))

    ## Packages in the db which do not satisfy the predicate.
    ind <- !predicate(db)
    ## Now find the recursive reverse dependencies of these and the
    ## non-standard packages missing from the db.
    rdepends <-
        tools:::package_dependencies(db1$Package[ind], db = db1,
                                     reverse = TRUE,
                                     recursive = recursive)
    rdepends <- unique(unlist(rdepends))
    ind[match(rdepends, db1$Package, nomatch = 0L)] <- TRUE

    ## And drop these from the db.
    db[!ind, , drop = FALSE]
}

available_packages_filters_db$`license/FOSS` <-
function(db) {
    predicate <- function(db)
        tools:::analyze_licenses(db[, "License"], db)$is_verified
    filter_packages_by_depends_predicates(db, predicate)
}

available_packages_filters_db$`license/restricts_use` <-
function(db) {
    predicate <- function(db) {
        ru <- tools:::analyze_licenses(db[, "License"], db)$restricts_use
        !is.na(ru) & !ru
    }
    filter_packages_by_depends_predicates(db, predicate)
}

available_packages_filters_db$CRAN <-
function(db)
{
    packages <- db[, "Package"]
    dups <- packages[duplicated(packages)]
    drop <- integer()
    CRAN <- getOption("repos")["CRAN"]
    for(d in dups) {
        pos <- which(packages == d)
        drop <- c(drop, pos[substring(db[pos, "Repository"], 1,
                                      nchar(CRAN)) != CRAN])
    }
    if(length(drop)) db[-drop, , drop = FALSE] else db
}


## unexported helper function
simplifyRepos <- function(repos, type)
{
    tail <- substring(contrib.url("---", type), 4L)
    ind <- regexpr(tail, repos, fixed=TRUE)
    ind <- ifelse(ind > 0L, ind-1L, nchar(repos, type="c"))
    substr(repos, 1L, ind)
}

update.packages <- function(lib.loc = NULL, repos = getOption("repos"),
                            contriburl = contrib.url(repos, type),
                            method, instlib = NULL, ask = TRUE,
                            available = NULL, oldPkgs = NULL, ...,
                            checkBuilt = FALSE, type = getOption("pkgType"))
{
    force(ask)  # just a check that it is valid before we start work
    text.select <- function(old)
    {
        update <- NULL
        for(k in seq_len(nrow(old))) {
            cat(old[k, "Package"], ":\n",
                "Version", old[k, "Installed"],
                "installed in", old[k, "LibPath"],
                if(checkBuilt) paste("built under R", old[k, "Built"]),
                "\n",
                "Version", old[k, "ReposVer"], "available at",
                simplifyRepos(old[k, "Repository"], type))
            cat("\n")
            answer <- substr(readline("Update (y/N/c)?  "), 1L, 1L)
            if(answer == "c" | answer == "C") {
                cat("cancelled by user\n")
                return(invisible())
            }
            if(answer == "y" | answer == "Y")
                update <- rbind(update, old[k,])
        }
        update
    }

    if(is.null(lib.loc))
        lib.loc <- .libPaths()

    if(is.null(available))
        available <- available.packages(contriburl = contriburl,
                                        method = method)

    if(!is.matrix(oldPkgs) && is.character(oldPkgs)) {
    	subset <- oldPkgs
    	oldPkgs <- NULL
    } else
    	subset <- NULL

    if(is.null(oldPkgs)) {
        ## since 'available' is supplied, 'contriburl' and 'method' are unused
	oldPkgs <- old.packages(lib.loc = lib.loc,
				contriburl = contriburl, method = method,
				available = available, checkBuilt = checkBuilt)
	## prune package versions which are invisible to require()
	if(!is.null(oldPkgs)) {
	    pkg <- 0L
	    while(pkg < nrow(oldPkgs)) {
		pkg <- pkg + 1L
		if(find.package(oldPkgs[pkg], lib.loc = lib.loc) !=
		   find.package(oldPkgs[pkg], lib.loc = oldPkgs[pkg,2])) {
		    warning(sprintf("package '%s' in library '%s' will not be updated",
				    oldPkgs[pkg], oldPkgs[pkg, 2]),
			    call. = FALSE, immediate. = TRUE)
		    oldPkgs <- oldPkgs[-pkg, , drop = FALSE]
		    pkg <- pkg - 1L
		}
	    }
	}
	if(is.null(oldPkgs))
	    return(invisible())
    } else if (!(is.matrix(oldPkgs) && is.character(oldPkgs)))
	stop("invalid 'oldPkgs'; must be a character vector or a result from old.packages()")

    if(!is.null(subset)) {
    	oldPkgs <- oldPkgs[ rownames(oldPkgs) %in% subset, ,drop=FALSE]
    	if (nrow(oldPkgs) == 0)
    	    return(invisible())
    }

    update <- if(is.character(ask) && ask == "graphics") {
        if(.Platform$OS.type == "windows" || .Platform$GUI == "AQUA"
           || (capabilities("tcltk") && capabilities("X11"))) {
            k <- select.list(oldPkgs[,1L], oldPkgs[,1L], multiple = TRUE,
                             title = "Packages to be updated", graphics = TRUE)
            oldPkgs[match(k, oldPkgs[,1L]), , drop=FALSE]
        } else text.select(oldPkgs)
    } else if(isTRUE(ask)) text.select(oldPkgs)
    else oldPkgs


    if(length(update)) {
        if(is.null(instlib)) instlib <-  update[, "LibPath"]
        ## do this a library at a time, to handle dependencies correctly.
        libs <- unique(instlib)
        for(l in libs)
            install.packages(update[instlib == l , "Package"], l,
                             contriburl = contriburl, method = method,
                             available = available, ..., type = type)
    }
}

old.packages <- function(lib.loc = NULL, repos = getOption("repos"),
                         contriburl = contrib.url(repos, type),
                         instPkgs = installed.packages(lib.loc = lib.loc),
                         method, available = NULL, checkBuilt = FALSE,
                         type = getOption("pkgType"))
{
    if(is.null(lib.loc))
        lib.loc <- .libPaths()
    if(!missing(instPkgs)) {
        ## actually we need rather more than this
        if(!is.matrix(instPkgs) || !is.character(instPkgs[, "Package"]))
            stop("ill-formed 'instPkgs' matrix")
    }
    if(NROW(instPkgs) == 0L) return(NULL)

    available <- if(is.null(available))
        available.packages(contriburl = contriburl, method = method)
    else tools:::.remove_stale_dups(available)

    update <- NULL

    currentR <- minorR <- getRversion()
    minorR[[c(1L, 3L)]] <- 0L # set patchlevel to 0
    for(k in 1L:nrow(instPkgs)) {
        if (instPkgs[k, "Priority"] %in% "base") next
        z <- match(instPkgs[k, "Package"], available[, "Package"])
        if(is.na(z)) next
        onRepos <- available[z, ]
        ## works OK if Built: is missing (which it should not be)
	if((!checkBuilt || package_version(instPkgs[k, "Built"]) >= minorR) &&
           package_version(onRepos["Version"]) <=
           package_version(instPkgs[k, "Version"])) next
        deps <- onRepos["Depends"]
        if(!is.na(deps)) {
            Rdeps <- tools:::.split_dependencies(deps)[["R", exact=TRUE]]
            if(length(Rdeps) > 1L) {
                target <- Rdeps$version
                res <- do.call(Rdeps$op, list(currentR, target))
 ##               res <- eval(parse(text=paste("currentR", Rdeps$op, "target")))
                if(!res) next
            }
        }
        update <- rbind(update,
                        c(instPkgs[k, c("Package", "LibPath", "Version", "Built")],
                          onRepos["Version"], onRepos["Repository"]))
    }
    if(!is.null(update))
        colnames(update) <- c("Package", "LibPath", "Installed", "Built",
                              "ReposVer", "Repository")
    rownames(update) <- update[, "Package"]
    ## finally, remove any duplicate rows
    update[!duplicated(update), , drop = FALSE]
}

new.packages <- function(lib.loc = NULL, repos = getOption("repos"),
                         contriburl = contrib.url(repos, type),
                         instPkgs = installed.packages(lib.loc = lib.loc),
                         method, available = NULL, ask = FALSE,
                         ..., type = getOption("pkgType"))
{
    ask  # just a check that it is valid before we start work
    if(is.null(lib.loc)) lib.loc <- .libPaths()
    if(!is.matrix(instPkgs))
        stop(gettextf("no installed packages for (invalid?) 'lib.loc=%s'",
                      lib.loc), domain = NA)
    if(is.null(available))
        available <- available.packages(contriburl = contriburl,
                                        method = method)

    installed <- unique(instPkgs[, "Package"])

    poss <- sort(unique(available[ ,"Package"])) # sort in local locale
    res <- setdiff(poss, installed)

    update <- character()
    graphics <- FALSE
    if(is.character(ask) && ask == "graphics") {
        ask <- TRUE
        if(.Platform$OS.type == "windows" || .Platform$GUI == "AQUA"
           || (capabilities("tcltk") && capabilities("X11")))
            graphics <- TRUE
    }
    if(isTRUE(ask)) {
        if(length(res))
            update <- res[match(select.list(res, multiple = TRUE,
                                            title = "New packages to be installed",
                                            graphics = graphics)
                                , res)]
        else message("no new packages are available")
    }
    if(length(update)) {
        install.packages(update, lib = lib.loc[1L], contriburl = contriburl,
                         method = method, available = available,
                         type = type, ...)
        # Now check if they were installed and update 'res'
        dirs <- list.files(lib.loc[1L])
        updated <- update[update %in% dirs]
        res <- res[!res %in% updated]
    }
    res
}

.instPkgFields <- function(fields) {
    ## to be used in installed.packages() and similar
    requiredFields <-
        c(tools:::.get_standard_repository_db_fields(), "Built")
    if (is.null(fields))
	fields <- requiredFields
    else {
	stopifnot(is.character(fields))
	fields <- unique(c(requiredFields, fields))
    }
    ## Don't retain 'Package' and 'LibPath' fields as these are used to
    ## record name and path of installed packages.
    fields[! fields %in% c("Package", "LibPath")]
}


## Read packages' Description and aggregate 'fields' into a character matrix
## NB: this does not handle encodings, so only suitable for ASCII-only fields.
.readPkgDesc <- function(lib, fields, pkgs = list.files(lib))
{
    ## to be used in installed.packages() and similar
    ## As from 2.13.0 only look at metadata.
    ret <- matrix(NA_character_, length(pkgs), 2L+length(fields))
    for(i in seq_along(pkgs)) {
        pkgpath <- file.path(lib, pkgs[i])
        if(file.access(pkgpath, 5L)) next
        if (file.exists(file <- file.path(pkgpath, "Meta", "package.rds"))) {
            ## this is vulnerable to installs going on in parallel
            md <- try(readRDS(file))
            if(inherits(md, "try-error")) next
            desc <- md$DESCRIPTION[fields]
            if (!length(desc)) {
                warning(gettextf("metadata of %s is corrupt", sQuote(pkgpath)),
                        domain = NA)
                next
            }
            if("Built" %in% fields) {
                ## This should not be missing.
                if(is.null(md$Built$R)) {
                    warning(gettextf("metadata of %s is corrupt",
                                     sQuote(pkgpath)), domain = NA)
                    next
                }
                desc["Built"] <- as.character(md$Built$R)
            }
            ret[i, ] <- c(pkgs[i], lib, desc)
        }
    }
    ret[!is.na(ret[, 1L]), ]
}

installed.packages <-
    function(lib.loc = NULL, priority = NULL, noCache = FALSE,
             fields = NULL, subarch = .Platform$r_arch)
{
    if(is.null(lib.loc))
        lib.loc <- .libPaths()
    if(!is.null(priority)) {
        if(!is.character(priority))
            stop("'priority' must be character or NULL")
        if(any(b <- priority %in% "high"))
            priority <- c(priority[!b], "recommended","base")
    }

    fields <- .instPkgFields(fields)
    retval <- matrix(character(), 0L, 2L + length(fields))
    for(lib in lib.loc) {
        if(noCache) {
            ret0 <- .readPkgDesc(lib, fields)
            if(length(ret0)) retval <- rbind(retval, ret0)
        } else {
            ## Previously used URLencode for e.g. Windows paths with drives
            ## This version works for very long file names.
            base <- paste(c(lib, fields), collapse = ",")
            ## add length and 64-bit CRC in hex (in theory, seems
            ## it is actually 32-bit on some systems)
            enc <- sprintf("%d_%s", nchar(base), .Call(C_crc64, base))
            dest <- file.path(tempdir(), paste0("libloc_", enc, ".rds"))
            if(file.exists(dest) &&
               file.info(dest)$mtime > file.info(lib)$mtime &&
               (val <- readRDS(dest))$base == base)
                ## use the cache file
                retval <- rbind(retval, val$value)
            else {
                ret0 <- .readPkgDesc(lib, fields)
                if(length(ret0)) {
                    retval <- rbind(retval, ret0)
                    ## save the cache file
                    saveRDS(list(base = base, value = ret0), dest)
                }
            }
        }
    }

    .fixupPkgMat(retval, fields, priority, subarch)
}

.fixupPkgMat <- function(mat, fields, priority, subarch=NULL)
{
    ## to be used in installed.packages() and similar
    colnames(mat) <- c("Package", "LibPath", fields)
    if (length(mat) && !is.null(priority)) {
	keep <- !is.na(pmatch(mat[,"Priority"], priority,
			      duplicates.ok = TRUE))
	mat <- mat[keep, , drop = FALSE]
    }
    if (length(mat) && !is.null(subarch) && nzchar(subarch)) {
        archs <- strsplit(mat[, "Archs"], ", ", fixed = TRUE)
        keep <- unlist(lapply(archs,
                              function(x) is.na(x[1L]) || subarch %in% x))
	mat <- mat[keep, , drop = FALSE]
    }
    if (length(mat)) mat <- mat[, colnames(mat) != "Archs", drop = FALSE]
    if (length(mat)) rownames(mat) <- mat[, "Package"]
    mat
}


remove.packages <- function(pkgs, lib)
{
    updateIndices <- function(lib) {
        ## This matches what install.packages() does
        if(lib == .Library && .Platform$OS.type == "unix") {
            message("Updating HTML index of packages in '.Library'")
            make.packages.html(.Library)
        }
        ## FIXME: only needed for packages installed < 2.13.0,
        ## so remove eventually
        ## is this the lib now empty?
        Rcss <- file.path(lib, "R.css")
        if (file.exists(Rcss)) {
            pkgs <- Sys.glob(file.path(lib, "*", "Meta", "package.rds"))
            if (!length(pkgs)) unlink(Rcss)
        }
    }

    if(!length(pkgs)) return(invisible())

    if(missing(lib) || is.null(lib)) {
        lib <- .libPaths()[1L]
	message(sprintf(ngettext(length(pkgs),
                                 "Removing package from %s\n(as %s is unspecified)",
                                 "Removing packages from %s\n(as %s is unspecified)"),
                        sQuote(lib), sQuote("lib")), domain = NA)
    }

    paths <- find.package(pkgs, lib)
    if(length(paths)) {
        unlink(paths, TRUE)
        for(lib in unique(dirname(paths))) updateIndices(lib)
    }
    invisible()
}

download.packages <- function(pkgs, destdir, available = NULL,
                              repos = getOption("repos"),
                              contriburl = contrib.url(repos, type),
                              method, type = getOption("pkgType"), ...)
{
    dirTest <- function(x) !is.na(isdir <- file.info(x)$isdir) & isdir

    nonlocalcran <- length(grep("^file:", contriburl)) < length(contriburl)
    if(nonlocalcran && !dirTest(destdir))
        stop("'destdir' is not a directory")
    if(is.null(available))
        available <- available.packages(contriburl=contriburl, method=method)

    retval <- matrix(character(), 0L, 2L)
    for(p in unique(pkgs))
    {
        ok <- (available[,"Package"] == p)
        ok <- ok & !is.na(ok)
        if(!any(ok))
            warning(gettextf("no package %s at the repositories", sQuote(p)),
                    domain = NA, immediate. = TRUE)
        else {
            if(sum(ok) > 1L) { # have multiple copies
                vers <- package_version(available[ok, "Version"])
                keep <- vers == max(vers)
                keep[duplicated(keep)] <- FALSE
                ok[ok][!keep] <- FALSE
            }
            if (substr(type, 1L, 10L) == "mac.binary") type <- "mac.binary"
            ## in Oct 2009 we introduced file names in PACKAGES files
            File <- available[ok, "File"]
            fn <- paste0(p, "_", available[ok, "Version"],
                         switch(type,
                                "source" = ".tar.gz",
                                "mac.binary" = ".tgz",
                                "win.binary" = ".zip"))
            have_fn <- !is.na(File)
            fn[have_fn] <- File[have_fn]
            repos <- available[ok, "Repository"]
            if(length(grep("^file:", repos)) > 0L) { # local repository
                ## This could be file: + file path or a file:/// URL.
                if(substring(repos, 1L, 8L) == "file:///") {
                    ## We need to derive the file name from the URL
                    ## This is tricky as so many forms have been allowed,
                    ## and indeed external methods may do even more.
                    fn <- paste(substring(repos, 8L), fn, sep = "/")
                    ## This leaves a path beginning with /
                    if(.Platform$OS.type == "windows") {
                        if(length(grep("^/[A-Za-z]:", fn)))
                            fn <- substring(fn, 2L)
                    }
                } else {
                    fn <- paste(substring(repos, 6L), fn, sep = "/")
                }
                if(file.exists(fn))
                    retval <- rbind(retval, c(p, fn))
                else
                    warning(gettextf("package %s does not exist on the local repository", sQuote(p)),
                            domain = NA, immediate. = TRUE)
            } else {
                url <- paste(repos, fn, sep = "/")
                destfile <- file.path(destdir, fn)

                res <- try(download.file(url, destfile, method, mode="wb", ...))
                if(!inherits(res, "try-error") && res == 0L)
                    retval <- rbind(retval, c(p, destfile))
                else
                    warning(gettextf("download of package %s failed", sQuote(p)),
                            domain = NA, immediate. = TRUE)
            }
        }
    }

    retval
}

contrib.url <- function(repos, type = getOption("pkgType"))
{
    ## Not entirely clear this is optimal
    if(type == "both") type <- "source"
    if(is.null(repos)) return(NULL)
    if("@CRAN@" %in% repos && interactive()) {
        cat(gettext("--- Please select a CRAN mirror for use in this session ---"),
            "\n", sep = "")
        flush.console()
        chooseCRANmirror()
        m <- match("@CRAN@", repos)
        nm <- names(repos)
        repos[m] <- getOption("repos")["CRAN"]
        if(is.null(nm)) nm <- rep("", length(repos))
        nm[m] <- "CRAN"
        names(repos) <- nm
    }
    if("@CRAN@" %in% repos) stop("trying to use CRAN without setting a mirror")

    ver <- paste(R.version$major,
                 strsplit(R.version$minor, ".", fixed=TRUE)[[1L]][1L], sep = ".")
    mac.path <- "macosx"
    if (substr(type, 1L, 11L) == "mac.binary.") {
        mac.path <- paste(mac.path, substring(type, 12L), sep = "/")
        type <- "mac.binary"
    }
    res <- switch(type,
		"source" = paste(gsub("/$", "", repos), "src", "contrib", sep = "/"),
                "mac.binary" = paste(gsub("/$", "", repos), "bin", mac.path, "contrib", ver, sep = "/"),
                "win.binary" = paste(gsub("/$", "", repos), "bin", "windows", "contrib", ver, sep = "/")
               )
    res
}


getCRANmirrors <- function(all = FALSE, local.only = FALSE)
{
    m <- NULL
    if(!local.only) {
        ## try to handle explicitly failure to connect to CRAN.
        con <- url("http://cran.r-project.org/CRAN_mirrors.csv")
        m <- try(open(con, "r"), silent = TRUE)
        if(!inherits(m, "try-error")) m <- try(read.csv(con, as.is = TRUE))
        close(con)
    }
    if(is.null(m) || inherits(m, "try-error"))
        m <- read.csv(file.path(R.home("doc"), "CRAN_mirrors.csv"),
                      as.is = TRUE)
    if(!all) m <- m[as.logical(m$OK), ]
    m
}


chooseCRANmirror <- function(graphics = getOption("menu.graphics"), ind = NULL)
{
    if(is.null(ind) && !interactive())
        stop("cannot choose a CRAN mirror non-interactively")
    m <- getCRANmirrors(all = FALSE, local.only = FALSE)
    res <- if (length(ind)) as.integer(ind)[1L] else
    menu(m[, 1L], graphics, "CRAN mirror")
    if(res > 0L) {
        URL <- m[res, "URL"]
        repos <- getOption("repos")
        repos["CRAN"] <- gsub("/$", "", URL[1L])
        options(repos = repos)
    }
    invisible()
}

chooseBioCmirror <- function(graphics = getOption("menu.graphics"), ind = NULL)
{
    if(is.null(ind) && !interactive())
        stop("cannot choose a BioC mirror non-interactively")
    m <- c("Seattle (USA)"="http://www.bioconductor.org"
	   , "Bethesda (USA)"="http://watson.nci.nih.gov/bioc_mirror"
	   , "Dortmund (Germany)"="http://bioconductor.statistik.tu-dortmund.de"
	   , "Anhui (China)"="http://mirrors.ustc.edu.cn/bioc/"
	   , "Cambridge (UK)"="http://mirrors.ebi.ac.uk/bioconductor/"
	   , "Riken, Kobe (Japan)" = "http://bioconductor.jp/"
	   , "Canberra (Australia)" = "http://mirror.aarnet.edu.au/pub/bioconductor/"
	   , "Sao Paulo (Brazil)" = "http://bioconductor.fmrp.usp.br/"
	   )
    res <- if (length(ind)) as.integer(ind)[1L] else
    menu(names(m), graphics, "BioC mirror")
    if(res > 0L) options("BioC_mirror" = m[res])
    invisible()
}

setRepositories <-
    function(graphics = getOption("menu.graphics"), ind = NULL,
             addURLs = character())
{
    if(is.null(ind) && !interactive())
        stop("cannot set repositories non-interactively")
    p <- file.path(Sys.getenv("HOME"), ".R", "repositories")
    if(!file.exists(p))
        p <- file.path(R.home("etc"), "repositories")
    a <- tools:::.read_repositories(p)
    pkgType <- getOption("pkgType")
    if (pkgType == "both") pkgType <- .Platform$pkgType
    if(length(grep("^mac\\.binary", pkgType))) pkgType <- "mac.binary"
    thisType <- a[[pkgType]]
    a <- a[thisType, 1L:3L]
    repos <- getOption("repos")
    ## Now look for CRAN and any others in getOptions("repos")
    if("CRAN" %in% row.names(a) && !is.na(CRAN <- repos["CRAN"]))
        a["CRAN", "URL"] <- CRAN
    ## Set as default any already in the option.
    a[(a[["URL"]] %in% repos), "default"] <- TRUE
    new <- !(repos %in% a[["URL"]])
    if(any(new)) {
        aa <- names(repos[new])
        if(is.null(aa)) aa <- rep("", length(repos[new]))
        aa[aa == ""] <- repos[new][aa == ""]
        newa <- data.frame(menu_name=aa, URL=repos[new], default=TRUE)
        row.names(newa) <- aa
        a <- rbind(a, newa)
    }

    default <- a[["default"]]

    res <- if(length(ind)) as.integer(ind)
    else {
        title <- if(graphics) "Repositories" else gettext("--- Please select repositories for use in this session ---\n")
        match(select.list(a[, 1L], a[default, 1L], multiple = TRUE, title,
                           graphics = graphics), a[, 1L])
    }
    if(length(res) || length(addURLs)) {
        repos <- a[["URL"]]
        names(repos) <- row.names(a)
        repos <- c(repos[res], addURLs)
        options(repos = repos)
    }
}



## used in some BioC packages and their support in tools.
compareVersion <- function(a, b)
{
    if(is.na(a)) return(-1L)
    if(is.na(b)) return(1L)
    a <- as.integer(strsplit(a, "[\\.-]")[[1L]])
    b <- as.integer(strsplit(b, "[\\.-]")[[1L]])
    for(k in seq_along(a))
        if(k <= length(b)) {
            if(a[k] > b[k]) return(1) else if(a[k] < b[k]) return(-1L)
        } else return(1L)
    if(length(b) > length(a)) return(-1L) else return(0L)
}

## ------------- private functions --------------------
.clean_up_dependencies <- function(x, available = NULL)
{
    ## x is a character vector of Depends / Suggests / Imports entries
    ## returns a character vector of all the package dependencies mentioned
    x <- x[!is.na(x)]
    if(!length(x)) return(x)
    x <- unlist(strsplit(x, ","))
    unique(sub("^[[:space:]]*([[:alnum:].]+).*$", "\\1" , x))
}

.clean_up_dependencies2 <- function(x, installed, available)
{
    ## x is a character vector of Depends / Suggests / Imports entries.
    ## Returns a list of length 2, a character vector of the names of
    ## all the package dependencies mentioned that are not already
    ## satisfied and one of those which cannot be satisfied (possibly
    ## of the form "pkg (>= ver)')

    .split_dependencies <- function(x) {
        .split2 <- function(x) {
            ## some have had space before ,
            x <- sub('[[:space:]]+$', '', x)
            x <- unique(sub("^[[:space:]]*(.*)", "\\1" , x))
            names(x) <- sub("^([[:alnum:].]+).*$", "\\1" , x)
            x <- x[names(x) != "R"]
	    x <- x[nzchar(x)]
            ## FIXME: a better way to handle duplicates.
            ## However, there should not be any, and if there are
            ## Depends: should be the first.
            x <- x[!duplicated(names(x))]
            lapply(x, tools:::.split_op_version)
        }
        ## given one of more concatenations of Depends/Imports/Suggests fields,
        ## return a named list of list(name, [op, version])
        if(!any(nzchar(x))) return(list())
        unlist(lapply(strsplit(x, ","), .split2), FALSE, FALSE)
    }
    x <- x[!is.na(x)]
    if(!length(x)) return(list(character(), character()))
    xx <- .split_dependencies(x)
    if(!length(xx)) return(list(character(), character()))
    ## Then check for those we already have installed
    pkgs <- installed[, "Package"]
    have <- sapply(xx, function(x) {
        if(length(x) == 3L) {
            if (! x[[1L]] %in% pkgs ) return(FALSE)
            if(x[[2L]] != ">=") return(TRUE)
            ## We may have the package installed more than once
            ## which we get will depend on the .libPaths() order,
            ## so for now just see if any installed version will do.
            current <- as.package_version(installed[pkgs == x[[1L]], "Version"])
            target <- as.package_version(x[[3L]])
            any(do.call(x$op, list(current, target)))
##            eval(parse(text = paste("any(current", x$op, "target)")))
        } else x[[1L]] %in% pkgs
    })
    xx <- xx[!have]
    if(!length(xx)) return(list(character(), character()))
    ## now check if we can satisfy the missing dependencies
    pkgs <- row.names(available)
    canget <- miss <- character()
    for (i in seq_along(xx)) {
        x <- xx[[i]]
        if(length(x) == 3L) {
            if (! x[[1L]] %in% pkgs ) { miss <- c(miss, x[[1L]]); next }
            if(x[[2L]] != ">=") { canget <- c(canget, x[[1L]]); next }
            ## we may have the package available more than once
            ## install.packages() will find the highest version.
            current <- as.package_version(available[pkgs == x[[1L]], "Version"])
            target <- as.package_version(x[[3L]])
            res <- any(do.call(x$op, list(current, target)))
##            res <- eval(parse(text = paste("any(current", x$op, "target)")))
            if(res) canget <- c(canget, x[[1L]])
            else  miss <- c(miss, paste0(x[[1L]], " (>= ", x[[3L]], ")"))
        } else if(x[[1L]] %in% pkgs) canget <- c(canget, x[[1L]])
        else miss <- c(miss, x[[1L]])
    }
    list(canget, miss)
}

.make_dependency_list <-
    function(pkgs, available,
             dependencies = c("Depends", "Imports", "LinkingTo"),
             recursive = FALSE)
{
    ## given a character vector of packages,
    ## return a named list of character vectors of their dependencies.
    ## If recursive = TRUE, do this recursively.
    if(!length(pkgs)) return(NULL)
    if(is.null(available))
        stop(gettextf("%s must be supplied", sQuote(available)), domain = NA)
    info <- available[pkgs, dependencies, drop = FALSE]
    x <- vector("list", length(pkgs)); names(x) <- pkgs
    if(recursive) {
        known <- row.names(available)
        xx <- vector("list", length(known)); names(xx) <- known
        info2 <-  available[, dependencies, drop = FALSE]
        for (i in seq_along(known))
            xx[[i]] <- .clean_up_dependencies(info2[i, ])
        for (i in pkgs) {
            p <- xx[[i]]
            p <- p[p %in% known]; p1 <- p
            repeat {
                extra <- unlist(xx[p1])
                extra <- extra[extra != i]
                extra <- extra[extra %in% known]
                deps <- unique(c(p, extra))
                if (length(deps) <= length(p)) break
                p1 <- deps[!deps %in% p]
                p <- deps
            }
            x[[i]] <- p
        }
    } else {
        for (i in seq_along(pkgs)) x[[i]] <- .clean_up_dependencies(info[i, ])
    }
    x
}

.find_install_order <- function(pkgs, dependencyList)
{
    ## given a character vector of packages, find an install order
    ## which reflects their dependencies.
    DL <- dependencyList[pkgs]
    ## some of the packages may be already installed, but the
    ## dependencies apply to those being got from CRAN.
    DL <- lapply(DL, function(x) x[x %in% pkgs])
    lens <- sapply(DL, length)
    if(all(lens > 0L)) {
        warning("every package depends on at least one other")
        return(pkgs)
    }
    done <- names(DL[lens == 0L]); DL <- DL[lens > 0L]
    while(length(DL)) {
        OK <- sapply(DL, function(x) all(x %in% done))
        if(!any(OK)) {
            warning(gettextf("packages %s are mutually dependent",
                             paste(sQuote(names(DL)), collapse = ", ")),
                    domain = NA)
            return(c(done,  names(DL)))
        }
        done <- c(done, names(DL[OK]))
        DL <- DL[!OK]
    }
    done
}
#  File src/library/utils/R/packages2.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

if (.Platform$OS.type == "windows")
    .install.macbinary <- function(...) NULL	# globalVariables isn't available, so use this to suppress the warning

getDependencies <-
    function(pkgs, dependencies = NA, available = NULL, lib = .libPaths()[1L])
{
    if (is.null(dependencies)) return(unique(pkgs))
    oneLib <- length(lib) == 1L
    dep2 <- NULL
    if(is.logical(dependencies) && is.na(dependencies))
        dependencies <- c("Depends", "Imports", "LinkingTo")
    depends <-
        is.character(dependencies) || (is.logical(dependencies) && dependencies)
    if(depends && is.logical(dependencies)) {
        dependencies <-  c("Depends", "Imports", "LinkingTo", "Suggests")
        dep2 <- c("Depends", "Imports", "LinkingTo")
    }
    if(depends && !oneLib) {
        warning("Do not know which element of 'lib' to install dependencies into\nskipping dependencies")
        depends <- FALSE
    }
    p0 <- unique(pkgs)
    miss <-  !p0 %in% row.names(available)
    if(sum(miss)) {
	warning(sprintf(ngettext(sum(miss),
				 "package %s is not available (for %s)",
				 "packages %s are not available (for %s)"),
			paste(sQuote(p0[miss]), collapse=", "),
			sub(" *\\(.*","", R.version.string)),
                domain = NA, call. = FALSE)
        if (sum(miss) == 1L &&
            !is.na(w <- match(tolower(p0[miss]),
                              tolower(row.names(available))))) {
            warning(sprintf("Perhaps you meant %s ?",
                            sQuote( row.names(available)[w])),
                    call. = FALSE, domain = NA)
        }
        flush.console()
    }
    p0 <- p0[!miss]

    if(depends) { # check for dependencies, recursively
        p1 <- p0 # this is ok, as 1 lib only
        ## INSTALL prepends 'lib' to the libpath
        ## Here we are slightly more conservative
        libpath <- .libPaths()
        if(!lib %in% libpath) libpath <- c(lib, libpath)
        installed <- installed.packages(lib.loc = libpath,
                                        fields = c("Package", "Version"))
        not_avail <- character()
	repeat {
	    deps <- apply(available[p1, dependencies, drop = FALSE],
                          1L, function(x) paste(x[!is.na(x)], collapse=", "))
	    res <- .clean_up_dependencies2(deps, installed, available)
            not_avail <- c(not_avail, res[[2L]])
            deps <- unique(res[[1L]])
            ## R should not get to here, but be safe
            deps <- deps[!deps %in% c("R", pkgs)]
	    if(!length(deps)) break
	    pkgs <- c(deps, pkgs)
	    p1 <- deps
            if(!is.null(dep2)) { dependencies <- dep2; dep2 <- NULL }
	}
        if(length(not_avail)) {
            not_avail <- unique(not_avail)
            warning(sprintf(ngettext(length(not_avail),
                                     "dependency %s is not available",
                                     "dependencies %s are not available"),
                            paste(sQuote(not_avail), collapse=", ")),
                    domain = NA, call. = FALSE, immediate. = TRUE)
            flush.console()
        }

        pkgs <- unique(pkgs)
        pkgs <- pkgs[pkgs %in% row.names(available)]
        if(length(pkgs) > length(p0)) {
            added <- setdiff(pkgs, p0)
            message(sprintf(ngettext(length(added),
                                     "also installing the dependency %s",
                                     "also installing the dependencies %s"),
                            paste(sQuote(added), collapse=", ")),
                    "\n", domain = NA)
            flush.console()
        }
        p0 <- pkgs
    }
    p0
}

install.packages <-
    function(pkgs, lib, repos = getOption("repos"),
             contriburl = contrib.url(repos, type),
             method, available = NULL, destdir = NULL, dependencies = NA,
             type = getOption("pkgType"),
             configure.args = getOption("configure.args"),
             configure.vars = getOption("configure.vars"),
             clean = FALSE, Ncpus = getOption("Ncpus", 1L),
	     verbose = getOption("verbose"),
             libs_only = FALSE, INSTALL_opts, quiet = FALSE,
             keep_outputs = FALSE,
             ...)
{
    if (is.logical(clean) && clean)
        clean <- "--clean"
    if(is.logical(dependencies) && is.na(dependencies))
        dependencies <- if(!missing(lib) && length(lib) > 1L) FALSE
        else c("Depends", "Imports", "LinkingTo")

    ## Compute the configuration arguments for a given package.
    ## If configure.args is an unnamed character vector, use that.
    ## If it is named, match the pkg name to the names of the character
    ## vector and if we get a match, use that element.
    ## Similarly, configure.args is a list(), match pkg to the names pkg
    ## and use that element, collapsing it into a single string.

    get_package_name <- function(pkg) {
        ## Since the pkg argument can be the name of a file rather than
        ## a regular package name, we have to clean that up.
        gsub("_\\.(zip|tar\\.gz)", "",
             gsub(.standard_regexps()$valid_package_version, "",
                  basename(pkg)))
    }

    getConfigureArgs <- function(pkg)
    {
        if(.Platform$OS.type == "windows") return(character())

        pkg <- get_package_name(pkg)

        if(length(pkgs) == 1L && length(configure.args) &&
           length(names(configure.args)) == 0L)
            return(paste0("--configure-args=",
                          shQuote(paste(configure.args, collapse = " "))))

        if (length(configure.args) && length(names(configure.args))
              && pkg %in% names(configure.args))
            config <- paste0("--configure-args=",
                             shQuote(paste(configure.args[[ pkg ]], collapse = " ")))
        else
            config <- character()

        config
    }

    getConfigureVars <- function(pkg)
    {
        if(.Platform$OS.type == "windows") return(character())

        pkg <- get_package_name(pkg)

        if(length(pkgs) == 1L && length(configure.vars) &&
           length(names(configure.vars)) == 0L)
            return(paste0("--configure-vars=",
                          shQuote(paste(configure.vars, collapse = " "))))

        if (length(configure.vars) && length(names(configure.vars))
              && pkg %in% names(configure.vars))
            config <- paste0("--configure-vars=",
                             shQuote(paste(configure.vars[[ pkg ]], collapse = " ")))
        else
            config <- character()

        config
    }

    get_install_opts <- function(pkg) {
        if(!length(INSTALL_opts))
            character()
        else
            paste(INSTALL_opts[[get_package_name(pkg)]], collapse = " ")
    }

    if(missing(pkgs) || !length(pkgs)) {
        if(!interactive()) stop("no packages were specified")
        ## if no packages were specified, use a menu
	if(.Platform$OS.type == "windows" || .Platform$GUI == "AQUA"
           || (capabilities("tcltk")
               && capabilities("X11") && suppressWarnings(tcltk:::.TkUp)) ) {
            ## this is the condition for a graphical select.list()
	} else
	    stop("no packages were specified")

        ## This will only offer the specified type.
	if(is.null(available))
	    available <- available.packages(contriburl = contriburl,
					    method = method)
	if(NROW(available)) {
            ## avoid duplicate entries in menus, since the latest available
            ## will be picked up
            ## sort in the locale, as R <= 2.10.1 did so
	    pkgs <- select.list(sort(unique(rownames(available))),
                                multiple = TRUE,
                                title = "Packages", graphics = TRUE)
	}
	if(!length(pkgs)) stop("no packages were specified")
    }

    if(missing(lib) || is.null(lib)) {
        lib <- .libPaths()[1L]
	if(!quiet && length(.libPaths()) > 1L)
	    message(sprintf(ngettext(length(pkgs),
                                     "Installing package into %s\n(as %s is unspecified)",
                                     "Installing packages into %s\n(as %s is unspecified)"),
                            sQuote(lib), sQuote("lib")), domain = NA)
    }

    ## check for writability by user
    ok <- file.info(lib)$isdir & (file.access(lib, 2) == 0)
    if(length(lib) > 1 && any(!ok))
        stop(sprintf(ngettext(sum(!ok),
                              "'lib' element %s is not a writable directory",
                              "'lib' elements %s are not writable directories"),
                     paste(sQuote(lib[!ok]), collapse=", ")), domain = NA)
    if(length(lib) == 1L && .Platform$OS.type == "windows") {
        ## file.access is unreliable on Windows, especially >= Vista.
        ## the only known reliable way is to try it
        ok <- file.info(lib)$isdir %in% TRUE # dir might not exist, PR#14311
        if(ok) {
            fn <- file.path(lib, paste("_test_dir", Sys.getpid(), sep = "_"))
            unlink(fn, recursive = TRUE) # precaution
            res <- try(dir.create(fn, showWarnings = FALSE))
            if(inherits(res, "try-error") || !res) ok <- FALSE
            else unlink(fn, recursive = TRUE)
        }
    }
    if(length(lib) == 1L && !ok) {
        warning(gettextf("'lib = \"%s\"' is not writable", lib),
                domain = NA, immediate. = TRUE)
        userdir <- unlist(strsplit(Sys.getenv("R_LIBS_USER"),
                                   .Platform$path.sep))[1L]
	if(interactive()) {
	    ask.yes.no <- function(msg) {
                ##' returns "no" for "no",  otherwise 'ans', a string
		msg <- gettext(msg)
		if(.Platform$OS.type == "windows") {
		    ans <- winDialog("yesno", sprintf(msg, sQuote(userdir)))
		    if(ans != "YES") "no" else ans
		} else {
		    ans <- readline(paste(sprintf(msg, userdir), " (y/n) "))
		    if(substr(ans, 1L, 1L) == "n") "no" else ans
		}
	    }
	    ans <- ask.yes.no("Would you like to use a personal library instead?")
	    if(identical(ans, "no")) stop("unable to install packages")

	    lib <- userdir
	    if(!file.exists(userdir)) {
		ans <- ask.yes.no("Would you like to create a personal library\n%s\nto install packages into?")
		if(identical(ans, "no")) stop("unable to install packages")
		if(!dir.create(userdir, recursive = TRUE))
                    stop(gettextf("unable to create %s", sQuote(userdir)),
                         domain = NA)
		.libPaths(c(userdir, .libPaths()))
	    }
	} else stop("unable to install packages")
    }

    lib <- normalizePath(lib)

    ## check if we should infer repos=NULL
    if(length(pkgs) == 1L && missing(repos) && missing(contriburl)) {
        if((type == "source" && length(grep("\\.tar.gz$", pkgs))) ||
           (type %in% "win.binary" && length(grep("\\.zip$", pkgs))) ||
           (substr(type, 1L, 10L) == "mac.binary"
            && length(grep("\\.tgz$", pkgs)))) {
            repos <- NULL
            message("inferring 'repos = NULL' from the file name")
        }
    }

# for testing .Platform$pkgType <- "mac.binary.leopard"
    ## Look at type == "both"
    if (type == "both") {
        ## NB it is only safe to use binary packages with a Mac OS X
        ## build that uses the same R foundation layout as CRAN since
        ## paths in DSOs are hard-coded.
        type2 <- .Platform$pkgType
        if (type2 == "source")
            stop("type == \"both\" can only be used on Windows or a CRAN build for Mac OS X")
        if(!missing(contriburl) || !is.null(available))
            stop("type == \"both\" cannot be used if 'available' or 'contriburl' is specified")
        if(is.null(repos))
            stop("type == \"both\" cannot be used with 'repos = NULL'")
        type <- "source"
        contriburl <- contrib.url(repos, "source")
        # The line above may have changed the repos option, so..
        if (missing(repos)) repos <- getOption("repos")
        available <-
            available.packages(contriburl = contriburl, method = method,
                               fields = "NeedsCompilation")
        pkgs <- getDependencies(pkgs, dependencies, available, lib)
        ## Now see what we can get as binary packages.
        av2 <- available.packages(contriburl = contrib.url(repos, type2),
                                  method = method)
        bins <- row.names(av2)
        bins <- pkgs[pkgs %in% bins]
        srcOnly <- pkgs[! pkgs %in% bins]
        binvers <- av2[bins, "Version"]
        hasSrc <-  !is.na(av2[bins, "Archs"])

        srcvers <- available[bins, "Version"]
        later <- as.numeric_version(binvers) < srcvers
        if(any(later)) {
            msg <- ngettext(sum(later),
                            "There is a binary version available but the source version is later",
                            "There are binary versions available but the source versions are later")
            cat("\n",
                paste(strwrap(msg, indent = 2, exdent = 2), collapse = "\n"),
                ":\n", sep = "")
            out <- data.frame(`binary` = binvers, `source` = srcvers,
                              `needs_compilation` =  hasSrc,
                              row.names = bins,
                              check.names = FALSE)[later, ]
            print(out)
            cat("\n")
            if(interactive() && any(later & hasSrc)) {
                msg <-
                    ngettext(sum(later & hasSrc),
                             "Do you want to install from sources the package which need compilation?",
                             "Do you want to install from sources the packages which need compilation?")
                message(msg, domain = NA)
                res <- readline("y/n: ")
                if(res != "y") later <- later & !hasSrc
            }
        }
        bins <- bins[!later]

        if(interactive() && length(srcOnly)) {
            nc <- !( available[srcOnly, "NeedsCompilation"] %in% "no" )
            s2 <- srcOnly[nc]
            if(length(s2)) {
                msg <-
                    ngettext(length(s2),
                             "Package which are only available in source form, and may need compilation of C/C++/Fortran",
                             "Packages which are only available in source form, and may need compilation of C/C++/Fortran")
                msg <- c(paste0(msg, ": "), sQuote(s2))
                msg <- strwrap(paste(msg, collapse = " "), exdent = 2)
                message(paste(msg, collapse = "\n"), domain = NA)
                message("Do you want to attempt to install these from sources?")
                res <- readline("y/n: ")
                if(res != "y") pkgs <- setdiff(pkgs, s2)
            }
        }

        if(length(bins)) {
            if(type2 == "win.binary")
                .install.winbinary(pkgs = bins, lib = lib,
                                   contriburl = contrib.url(repos, type2),
                                   method = method, available = av2,
                                   destdir = destdir,
                                   dependencies = NULL,
                                   libs_only = libs_only, ...)
            else
                .install.macbinary(pkgs = bins, lib = lib,
                                   contriburl = contrib.url(repos, type2),
                                   method = method, available = av2,
                                   destdir = destdir,
                                   dependencies = NULL, ...)
        }
        pkgs <- setdiff(pkgs, bins)
        if(!length(pkgs)) return(invisible())
        message(sprintf(ngettext(length(pkgs),
                                     "installing the source package %s",
                                     "installing the source packages %s"),
                        paste(sQuote(pkgs), collapse=", ")),
                "\n", domain = NA)
	flush.console()
    } else if (getOption("install.packages.check.source", "yes") %in% "yes"
               && (type %in% "win.binary" || substr(type, 1L, 10L) == "mac.binary")) {
        if (missing(contriburl) && is.null(available) && !is.null(repos)) {
            contriburl2 <- contrib.url(repos, "source")
	    # The line above may have changed the repos option, so..
            if (missing(repos)) repos <- getOption("repos")
            av1 <- try(suppressWarnings(available.packages(contriburl = contriburl2, method = method)), silent = TRUE)
            if(inherits(av1, "try-error")) {
                message("source repository is unavailable to check versions")
                available <-
                    available.packages(contriburl = contrib.url(repos, type), method = method)
            } else {
                srcpkgs <- pkgs[pkgs %in% row.names(av1)]
                ## Now see what we can get as binary packages.
                available <-
                    available.packages(contriburl = contrib.url(repos, type), method = method)
                bins <- pkgs[pkgs %in% row.names(available)]
                ## so a package might only be available as source,
                ## or it might be later in source.
                ## FIXME: might only want to check on the same repository,
                ## allowing for CRANextras.
                na <- srcpkgs[!srcpkgs %in% bins]
                if (length(na)) {
                    msg <-
                        sprintf(ngettext(length(na),
                                         "package %s is available as a source package but not as a binary",
                                         "packages %s are available as source packages but not as binaries"),
                                paste(sQuote(na), collapse = ", "))
                    cat("\n   ", msg, "\n\n", sep = "")
                }
                binvers <- available[bins, "Version"]
                srcvers <- binvers
                OK <- bins %in% srcpkgs
                srcvers[OK] <- av1[bins[OK], "Version"]
                later <- as.numeric_version(binvers) < srcvers
                if(any(later)) {
                    msg <- ngettext(sum(later),
                                    "There is a binary version available (and will be installed) but the source version is later",
                                    "There are binary versions available (and will be installed) but the source versions are later")
                    cat("\n",
                        paste(strwrap(msg, indent = 2, exdent = 2), collapse = "\n"),
                        ":\n", sep = "")
                    print(data.frame(`binary` = binvers, `source` = srcvers,
                                     row.names = bins,
                                     check.names = FALSE)[later, ])
                    cat("\n")
                }
            }
        }
    }

    if(.Platform$OS.type == "windows") {
        if(substr(type, 1L, 10L) == "mac.binary")
            stop("cannot install MacOS X binary packages on Windows")

        if(type %in% "win.binary") {
            ## include local .zip files
            .install.winbinary(pkgs = pkgs, lib = lib, contriburl = contriburl,
                               method = method, available = available,
                               destdir = destdir,
                               dependencies = dependencies,
                               libs_only = libs_only, quiet = quiet,  ...)
            return(invisible())
        }
        ## Avoid problems with spaces in pathnames.
        have_spaces <- grep(" ", pkgs)
        if(length(have_spaces)) {
            ## we want the short name for the directory,
            ## but not for a .tar.gz, and package names never contain spaces.
            p <- pkgs[have_spaces]
            dirs <- shortPathName(dirname(p))
            pkgs[have_spaces] <- file.path(dirs, basename(p))
        }
        ## Avoid problems with backslashes
        ## -- will mess up UNC names, but they don't work
        pkgs <- gsub("\\\\", "/", pkgs)
    } else {
        if(substr(type, 1L, 10L) == "mac.binary") {
            if(!length(grep("darwin", R.version$platform)))
                stop("cannot install MacOS X binary packages on this platform")
            .install.macbinary(pkgs = pkgs, lib = lib, contriburl = contriburl,
                               method = method, available = available,
                               destdir = destdir,
                               dependencies = dependencies, quiet = quiet, ...)
            return(invisible())
        }

        if(type %in% "win.binary")
            stop("cannot install Windows binary packages on this platform")

        if(!file.exists(file.path(R.home("bin"), "INSTALL")))
            stop("This version of R is not set up to install source packages\nIf it was installed from an RPM, you may need the R-devel RPM")
    }

    ## we need to ensure that R CMD INSTALL runs with the same
    ## library trees as this session.
    ## FIXME: At least on Windows, either run sub-R directly (to avoid sh)
    ## or run the install in the current process.
    libpath <- .libPaths()
    libpath <- libpath[! libpath %in% .Library]
    if(length(libpath))
        libpath <- paste(libpath, collapse = .Platform$path.sep)

    cmd0 <- file.path(R.home("bin"), "R")
    args0 <- c("CMD", "INSTALL")

    output <- if(quiet) FALSE else ""
    env <- character()

    outdir <- getwd()
    if(is.logical(keep_outputs)) {
        if(is.na(keep_outputs))
            keep_outputs <- FALSE
    } else if(is.character(keep_outputs) &&
              (length(keep_outputs) == 1L)) {
        if(!file_test("-d", keep_outputs) &&
           !dir.create(keep_outputs, recursive = TRUE))
            stop(gettextf("unable to create %s", sQuote(keep_outputs)),
                 domain = NA)
        outdir <- normalizePath(keep_outputs)
        keep_outputs <- TRUE
    } else
        stop(gettextf("invalid %s argument", sQuote("keep_outputs")),
             domain = NA)

    if(length(libpath)) {
        ## <NOTE>
        ## For the foreseeable future, the 'env' argument to system2()
        ## on Windows is limited to calls to make and rterm (but not R
        ## CMD): hence need to set the R_LIBS env var here.
        if(.Platform$OS.type == "windows") {
            ## We don't have a way to set an environment variable for
            ## a single command, as we do not spawn a shell.
            oldrlibs <- Sys.getenv("R_LIBS")
            Sys.setenv(R_LIBS = libpath)
            on.exit(Sys.setenv(R_LIBS = oldrlibs))
        } else
            env <- paste("R_LIBS", shQuote(libpath), sep = "=")
        ## </NOTE>
    }

    if (is.character(clean))
        args0 <- c(args0, clean)
    if (libs_only)
        args0 <- c(args0, "--libs-only")
    if (!missing(INSTALL_opts)) {
        if(!is.list(INSTALL_opts)) {
            args0 <- c(args0, paste(INSTALL_opts, collapse = " "))
            INSTALL_opts <- list()
        }
    } else {
        INSTALL_opts <- list()
    }

    if(verbose)
        message(gettextf("system (cmd0): %s",
                         paste(c(cmd0, args0), collapse = " ")),
                domain = NA)

    if(is.null(repos) & missing(contriburl)) {
        ## install from local source tarball(s)
        update <- cbind(path.expand(pkgs), lib) # for side-effect of recycling to same length

        for(i in seq_len(nrow(update))) {
            args <- c(args0,
                      get_install_opts(update[i, 1L]),
                      "-l", shQuote(update[i, 2L]),
                      getConfigureArgs(update[i, 1L]),
                      getConfigureVars(update[i, 1L]),
                      shQuote(update[i, 1L]))
            status <- system2(cmd0, args, env = env,
                              stdout = output, stderr = output)
            if(status > 0L)
                warning(gettextf("installation of package %s had non-zero exit status",
                                 sQuote(update[i, 1L])),
                        domain = NA)
            else if(verbose) {
                cmd <- paste(c(cmd0, args), collapse = " ")
                message(sprintf("%d): succeeded '%s'", i, cmd),
                        domain = NA)
            }
        }
        return(invisible())
    }

    tmpd <- destdir
    nonlocalrepos <- length(grep("^file:", contriburl)) < length(contriburl)
    if(is.null(destdir) && nonlocalrepos) {
        tmpd <- file.path(tempdir(), "downloaded_packages")
        if (!file.exists(tmpd) && !dir.create(tmpd))
            stop(gettextf("unable to create temporary directory %s",
                          sQuote(tmpd)),
                 domain = NA)
    }

    if(is.null(available))
        available <- available.packages(contriburl = contriburl,
                                        method = method)
    pkgs <- getDependencies(pkgs, dependencies, available, lib)

    foundpkgs <- download.packages(pkgs, destdir = tmpd, available = available,
                                   contriburl = contriburl, method = method,
                                   type = "source", quiet = quiet, ...)

    ## at this point 'pkgs' may contain duplicates,
    ## the same pkg in different libs
    if(length(foundpkgs)) {
	if(verbose) message(gettextf("foundpkgs: %s",
                                     paste(foundpkgs, collapse=", ")),
                            domain = NA)
        update <- unique(cbind(pkgs, lib))
        colnames(update) <- c("Package", "LibPath")
        found <- pkgs %in% foundpkgs[, 1L]
        files <- foundpkgs[match(pkgs[found], foundpkgs[, 1L]), 2L]
	if(verbose) message(gettextf("files: %s",
                                     paste(files, collapse=", \n\t")),
                            domain = NA)
        update <- cbind(update[found, , drop=FALSE], file = files)
        if(nrow(update) > 1L) {
            upkgs <- unique(pkgs <- update[, 1L])
            DL <- .make_dependency_list(upkgs, available)
            p0 <- .find_install_order(upkgs, DL)
            ## can't use update[p0, ] due to possible multiple matches
            update <- update[sort.list(match(pkgs, p0)), ]
        }

        if (Ncpus > 1L && nrow(update) > 1L) {
            ## if --no-lock or --lock was specified in INSTALL_opts
            ## that will override this.
            args0 <- c(args0, "--pkglock")
            tmpd <- file.path(tempdir(), "make_packages")
            if (!file.exists(tmpd) && !dir.create(tmpd))
                stop(gettextf("unable to create temporary directory %s",
                              sQuote(tmpd)),
                     domain = NA)
            mfile <- file.path(tmpd, "Makefile")
            conn <- file(mfile, "wt")
            deps <- paste(paste0(update[, 1L], ".ts"), collapse=" ")
            deps <- strwrap(deps, width = 75, exdent = 2)
            deps <- paste(deps, collapse=" \\\n")
            cat("all: ", deps, "\n", sep = "", file = conn)
            aDL <- .make_dependency_list(upkgs, available, recursive = TRUE)
            for(i in seq_len(nrow(update))) {
                pkg <- update[i, 1L]
                args <- c(args0,
                          get_install_opts(update[i, 3L]),
                          "-l", shQuote(update[i, 2L]),
                          getConfigureArgs(update[i, 3L]),
                          getConfigureVars(update[i, 3L]),
                          shQuote(update[i, 3L]),
                          ">", paste0(pkg, ".out"),
                          "2>&1")
                ## <NOTE>
                ## We currently only use env on Unix for R_LIBS.
                ## Windows we do Sys.setenv(R_LIBS = libpath),
                ## since system2() has limited support for 'env'
                ## Should we use env on Windows as well?
                ## If so, would we need
                ##   cmd <- paste(c(shQuote(command), env, args),
                ##                collapse = " ")
                ## on Windows?
                cmd <- paste(c(shQuote(cmd0), args), collapse = " ")
                ## </NOTE>
                deps <- aDL[[pkg]]
                deps <- deps[deps %in% upkgs]
                ## very unlikely to be too long
                deps <- if(length(deps))
                    paste(paste0(deps, ".ts"), collapse = " ") else ""
                cat(paste0(pkg, ".ts: ", deps),
                    paste("\t@echo begin installing package", sQuote(pkg)),
                    paste0("\t@", cmd, " && touch ", pkg, ".ts"),
                    paste0("\t@cat ", pkg, ".out"),
                    "", sep = "\n", file = conn)
            }
            close(conn)
            cwd <- setwd(tmpd)
            on.exit(setwd(cwd))
            ## MAKE will be set by sourcing Renviron
            status <- system2(Sys.getenv("MAKE", "make"),
                              c("-k -j", Ncpus),
                              stdout = output, stderr = output,
                              env = env)
            if(status > 0L) {
                ## Try to figure out which
                pkgs <- update[, 1L]
                tss <- sub("\\.ts$", "", dir(".", pattern = "\\.ts$"))
                failed <- pkgs[!pkgs %in% tss]
		for (pkg in failed) system(paste0("cat ", pkg, ".out"))
                warning(gettextf("installation of one of more packages failed,\n  probably %s",
                                 paste(sQuote(failed), collapse = ", ")),
                        domain = NA)
            }
            if(keep_outputs)
                file.copy(paste0(update[, 1L], ".out"), outdir)
            setwd(cwd); on.exit()
            unlink(tmpd, recursive = TRUE)
        } else {
            for(i in seq_len(nrow(update))) {
                outfile <- if(keep_outputs) {
                    paste0(update[i, 1L], ".out")
                } else output
                args <- c(args0,
                          get_install_opts(update[i, 3L]),
                          "-l", shQuote(update[i, 2L]),
                          getConfigureArgs(update[i, 3L]),
                          getConfigureVars(update[i, 3L]),
                          update[i, 3L])
                status <- system2(cmd0, args, env = env,
                                  stdout = outfile, stderr = outfile)
                if(!quiet && keep_outputs)
                    writeLines(readLines(outfile))
                if(status > 0L)
                    warning(gettextf("installation of package %s had non-zero exit status",
                                     sQuote(update[i, 1L])),
                            domain = NA)
		else if(verbose) {
                    cmd <- paste(c(cmd0, args), collapse = " ")
                    message(sprintf("%d): succeeded '%s'", i, cmd),
                            domain = NA)
                }
            }
            if(keep_outputs && (outdir != getwd()))
                file.copy(paste0(update[, 1L], ".out"), outdir)
        }
        if(!quiet && nonlocalrepos && !is.null(tmpd) && is.null(destdir))
            cat("\n", gettextf("The downloaded source packages are in\n\t%s",
                               sQuote(normalizePath(tmpd, mustWork = FALSE))),
                "\n", sep = "")
        ## update packages.html on Unix only if .Library was installed into
        libs_used <- unique(update[, 2L])
        if(.Platform$OS.type == "unix" && .Library %in% libs_used) {
            message("Updating HTML index of packages in '.Library'")
            make.packages.html(.Library)
        }
    } else if(!is.null(tmpd) && is.null(destdir)) unlink(tmpd, TRUE)

    invisible()
}

## treat variables as global in a package, for codetools & check
globalVariables <- function(names, package, add = TRUE) {
    .listFile <- ".__global__"
    .simplePackageName <- function(env) {
        if(exists(".packageName", envir = env, inherits = FALSE))
           get(".packageName", envir = env)
        else
            "(unknown package)"
    }
    if(missing(package)) {
        env <- topenv(parent.frame())
        package <- .simplePackageName(env)
    }
    else if(is.environment(package)) {
        env <- package
        package <- .simplePackageName(env)
    }
    else
        env <- asNamespace(package)
    if(exists(.listFile, envir = env, inherits = FALSE))
        current <- get(.listFile, envir = env)
    else
        current <- character()
    if(! missing(names)) {
        if(environmentIsLocked(env))
            stop(gettextf("The namespace for package \"%s\" is locked; no changes in the global variables list may be made.",
                          package))
        if(add)
            current <- unique(c(current, names))
        else
            current <- names
        assign(.listFile, current, envir = env)
    }
    current
}

packageName <- function(env = parent.frame()) {
    if (!is.environment(env)) stop("'env' must be an environment")
    env <- topenv(env)
    if (exists(".packageName", envir = env, inherits = FALSE))
	get(".packageName", envir = env, inherits = FALSE)
    else if (identical(env, .BaseNamespaceEnv))
	"base"
    else
	NULL
}
#  File src/library/utils/R/page.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

page <- function(x, method = c("dput", "print"), ...)
{
    ## local functions to parcel out '...'
    local.file.show <- function(file, title = subx, delete.file = TRUE,
                                pager = getOption("pager"), ...)
        file.show(file, title = title, delete.file = delete.file, pager = pager)
    local.dput <- function(x, file, title, delete.file, pager, ...)
        dput(x, file, ...)
    local.print <- function(x, title, delete.file, pager, ...)
        print(x, ...)

    if(is.character(x) && length(x) == 1L) {
        subx <- x
        parent <- parent.frame()
        if(exists(subx, envir = parent)) # inherits=TRUE is default
            x <- get(subx, envir = parent)
        else
            stop(gettextf("no object named '%s' to show", x), domain = NA)
    } else {
        subx <- deparse(substitute(x))
    }
    file <- tempfile("Rpage.")
    if(match.arg(method) == "dput")
        local.dput(x, file, ...)
    else {
        sink(file)
        local.print(x, ...)
        sink()
    }
    local.file.show(file, ...)
}
#  File src/library/utils/R/progressBar.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

txtProgressBar <-
    function(min = 0, max = 1, initial = 0, char = "=",
             width = NA, title, label, style = 1, file = "")
{
    if(!identical(file, "") &&
       !(inherits(file, "connection") && isOpen(file)))
        stop("'file' must be \"\" or an open connection object")
    if(! style %in% 1L:3L) style <- 1
    .val <- initial
    .killed <- FALSE
    .nb <- 0L
    .pc <- -1L # This ensures the initial value is displayed for style = 3
    nw <- nchar(char, "w")
    if(is.na(width)) {
        width <- getOption("width")
        if(style == 3L) width <- width - 10L
        width <- trunc(width/nw)
    }
    if (max <= min) stop("must have 'max' > 'min'")

    up1 <- function(value) {
        if(!is.finite(value) || value < min || value > max) return()
        .val <<- value
        nb <- round(width*(value - min)/(max - min))
        if(.nb < nb) {
            cat(paste(rep.int(char, nb-.nb), collapse=""), file = file)
            flush.console()
        } else if (.nb > nb) {
            cat("\r", paste(rep.int(" ", .nb*nw), collapse=""),
                "\r", paste(rep.int(char, nb), collapse=""),
                sep = "", file = file)
            flush.console()
        }
        .nb <<- nb
    }

    up2 <- function(value) {
        if(!is.finite(value) || value < min || value > max) return()
        .val <<- value
        nb <- round(width*(value - min)/(max - min))
        if(.nb <= nb) {
            cat("\r", paste(rep.int(char, nb), collapse=""),
                sep = "", file = file)
            flush.console()
        } else {
            cat("\r", paste(rep.int(" ", .nb*nw), collapse=""),
                "\r", paste(rep.int(char, nb), collapse=""),
                sep = "", file = file)
            flush.console()
        }
        .nb <<- nb
    }

    up3 <- function(value) {
        if(!is.finite(value) || value < min || value > max) return()
        .val <<- value
        nb <- round(width*(value - min)/(max - min))
        pc <- round(100*(value - min)/(max - min))
        if(nb == .nb && pc == .pc) return()
        cat(paste(c("\r  |", rep.int(" ", nw*width+6)), collapse=""),
            file = file)
        cat(paste(c("\r  |",
                    rep.int(char, nb),
                    rep.int(" ", nw*(width-nb)),
                    sprintf("| %3d%%", pc)
                    ), collapse=""), file = file)
        flush.console()
        .nb <<- nb
        .pc <<- pc
    }

    getVal <- function() .val
    kill <- function()
        if(!.killed) {
            cat("\n", file = file)
            flush.console()
            .killed <<- TRUE
        }
    up <- switch(style, up1, up2, up3)
    up(initial) # will check if in range
    structure(list(getVal=getVal, up=up, kill=kill),
              class = "txtProgressBar")
}

getTxtProgressBar <- function(pb)
{
    if(!inherits(pb, "txtProgressBar"))
       stop(gettextf("'pb' is not from class %s",
                     dQuote("txtProgressBar")),
            domain = NA)
    pb$getVal()
}

setTxtProgressBar <- function(pb, value, title = NULL, label = NULL)
{
    if(!inherits(pb, "txtProgressBar"))
        stop(gettextf("'pb' is not from class %s",
                      dQuote("txtProgressBar")),
             domain = NA)
    oldval <- pb$getVal()
    pb$up(value)
    invisible(oldval)
}

close.txtProgressBar <- function(con, ...)
{
    con$kill()
    invisible(NULL)
}
#  File src/library/utils/R/prompt.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

prompt <-
function(object, filename = NULL, name = NULL, ...)
    UseMethod("prompt")

prompt.default <-
function(object, filename = NULL, name = NULL,
         force.function = FALSE, ...)
{
    is.missing.arg <- function(arg)
        typeof(arg) == "symbol" && deparse(arg) == ""

    if(missing(name))
        name <- if(is.character(object))
            object
        else {
            name <- substitute(object)
            ## <FIXME>
            ## This used to be:
            ##     if(is.language(name) && !is.name(name))
            ##         name <- eval(name)
            ##     as.character(name)
            ## but what is this trying to do?
            ## It seems that the eval() will typically give the given
            ## object, and surely we cannot use that as the name (even
            ## if the subsequent as.character() does not fail ...)
            ## Better to be defensive about this, and handle only cases
            ## we know will make sense ...
            if(is.name(name))
                as.character(name)
            else if(is.call(name)
                    && (as.character(name[[1L]]) %in%
                        c("::", ":::", "getAnywhere"))) {
                name <- as.character(name)
                name[length(name)]
            }
            else
                stop("cannot determine a usable name")
            ## </FIXME>
        }

    if(is.null(filename))
        filename <- paste0(name, ".Rd")

    x <- if(!missing(object))
        object
    else {
        ## Better than get(); works when called in fun :
        x <- get(name, envir = parent.frame())
    }

    ## <FIXME>
    ## If not a function or forced to document a function (?), always
    ## assume data set.
    if(!(is.function(x) || force.function))
        return(promptData(x, filename = filename, name = name))
    ## </FIXME>

    n <- length(argls <- formals(x))
    if(n) {
        arg.names <- arg.n <- names(argls)
        arg.n[arg.n == "..."] <- "\\dots"
    }
    ## Construct the 'call' for \usage.
    Call <- paste0(name, "(")
    for(i in seq_len(n)) {                       # i-th argument
        Call <- paste0(Call, arg.names[i],
                       if(!is.missing.arg(argls[[i]]))
                       paste0(" = ",
                              ## need to backtick symbols
                              paste(deparse(argls[[i]],
                                            backtick = TRUE,
                                            width.cutoff = 500L),
                                    collapse="\n")))
        if(i != n) Call <- paste0(Call, ", ")
    }

    ## Construct the definition for \examples.
    x.def <- attr(x, "source")
    if(is.null(x.def))
        x.def <- deparse(x)
    if(any(br <- substr(x.def, 1L, 1L) == "}"))
        x.def[br] <- paste(" ", x.def[br])

    ## escape "%" :
    x.def <- gsub("%", "\\\\%", x.def)

    Rdtxt <-
        list(name = paste0("\\name{", name, "}"),
#             version = "\\Rdversion{1.1}",
             aliases = c(paste0("\\alias{", name, "}"),
             paste("%- Also NEED an '\\alias' for EACH other topic",
                   "documented here.")),
             title = "\\title{\n%%  ~~function to do ... ~~\n}",
             description = c("\\description{",
             paste("%%  ~~ A concise (1-5 lines) description of what",
                   "the function does. ~~"),
             "}"),
             usage = c("\\usage{", paste0(Call, ")"), "}",
             paste("%- maybe also 'usage' for other objects",
                   "documented here.")),
             arguments = NULL,
             details = c("\\details{",
             paste("%%  ~~ If necessary, more details than the",
                   "description above ~~"),
             "}"),
             value = c("\\value{",
             "%%  ~Describe the value returned",
             "%%  If it is a LIST, use",
             "%%  \\item{comp1 }{Description of 'comp1'}",
             "%%  \\item{comp2 }{Description of 'comp2'}",
             "%% ...",
             "}"),
             references = paste("\\references{\n%% ~put references to the",
             "literature/web site here ~\n}"),
             author = "\\author{\n%%  ~~who you are~~\n}",
             note = c("\\note{\n%%  ~~further notes~~\n}",
             "",
             paste("%% ~Make other sections like Warning with",
                   "\\section{Warning }{....} ~"),
             ""),
             seealso = paste("\\seealso{\n%% ~~objects to See Also as",
             "\\code{\\link{help}}, ~~~\n}"),
             examples = c("\\examples{",
             "##---- Should be DIRECTLY executable !! ----",
             "##-- ==>  Define data, use random,",
             "##--	or do  help(data=index)  for the standard data sets.",
             "",
             "## The function is currently defined as",
             x.def,
             "}"),
             keywords = c(paste("% Add one or more standard keywords,",
             "see file 'KEYWORDS' in the"),
             "% R documentation directory.",
             "\\keyword{ ~kwd1 }",
             "\\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line"))

    Rdtxt$arguments <- if(n)
        c("\\arguments{",
          paste0("  \\item{", arg.n, "}{",
                 "\n%%     ~~Describe \\code{", arg.n, "} here~~\n}"),
          "}") ## else NULL

    if(is.na(filename)) return(Rdtxt)

    cat(unlist(Rdtxt), file = filename, sep = "\n")

    message(gettextf("Created file named %s.", sQuote(filename)),
            "\n",
            gettext("Edit the file and move it to the appropriate directory."),
            domain = NA)

    invisible(filename)
}

prompt.data.frame <-
function(object, filename = NULL, name = NULL, ...)
{
    if(missing(name))
        name <-
            if(is.character(object))
                object
            else {
                name <- substitute(object)
                if(is.name(name))
                    as.character(name)
                else
                    stop("cannot determine a usable name")
            }
    if(is.null(filename))
        filename <- paste0(name, ".Rd")

    x <- if(!missing(object))
        object
    else {
        ## Better than get(); works when called in fun :
        x <- get(name, envir = parent.frame())
    }

    ## <FIXME>
    ## Always assume data set ???
    promptData(x, filename = filename, name = name)
    ## </FIXME>
}

promptData <-
function(object, filename = NULL, name = NULL)
{
    if(missing(name))
        name <-
            if(is.character(object))
                object
            else {
                name <- substitute(object)
                if(is.name(name))
                    as.character(name)
                else
                    stop("cannot determine a usable name")
            }
    if(is.null(filename))
        filename <- paste0(name, ".Rd")

    x <- if(!missing(object))
        object
    else {
        ## Better than get(); works when called in fun :
        x <- get(name, envir = parent.frame())
    }

    ## Construct the format.
    if(is.data.frame(x)) {
        fmt <- c("\\format{",
                 paste("  A data frame with",
                       nrow(x),
                       "observations on the following",
                       ifelse(ncol(x) == 1,
                              "variable.",
                              paste(ncol(x), "variables."))),
                 "  \\describe{")
        for(i in names(x)) {
            xi <- x[[i]]
            fmt <-
                c(fmt,
                  paste0("    \\item{\\code{", i, "}}{",
                         if(inherits(xi, "ordered")) {
                             paste("an", data.class(xi),
                                   "factor with levels",
                                   paste0("\\code{", levels(xi), "}",
                                          collapse = " < "),
                                   collapse = " ")
                         } else if(inherits(xi, "factor")) {
                             paste("a factor with levels",
                                   paste0("\\code{", levels(xi), "}",
                                          collapse = " "),
                                   collapse = " ")
                         } else if(is.vector(xi)) {
                             paste("a", data.class(xi), "vector")
                         } else if(is.matrix(xi)) {
                             paste("a matrix with", ncol(xi), "columns")
                         } else {
                             paste("a", data.class(xi))
                         },
                         "}"))
        }
        fmt <- c(fmt, "  }", "}")
    }
    else {
        tf <- tempfile(); on.exit(unlink(tf))
        sink(tf) ; str(object) ; sink()
        fmt <- c("\\format{",
                 "  The format is:",
                 scan(tf, "", quiet = !getOption("verbose"), sep = "\n"),
                 "}")
    }

    Rdtxt <-
        list(name = paste0("\\name{", name, "}"),
#             version = "\\Rdversion{1.1}",
             aliases = paste0("\\alias{", name, "}"),
             docType = "\\docType{data}",
             title = "\\title{\n%%   ~~ data name/kind ... ~~\n}",
             description = c("\\description{",
             "%%  ~~ A concise (1-5 lines) description of the dataset. ~~",
             "}"),
             usage = paste0("\\usage{data(", name, ")}"),
             format = fmt,
             details = c("\\details{",
             paste("%%  ~~ If necessary, more details than the",
                   "__description__ above ~~"),
             "}"),
             source = c("\\source{",
             paste("%%  ~~ reference to a publication or URL",
                   "from which the data were obtained ~~"),
             "}"),
             references = c("\\references{",
             "%%  ~~ possibly secondary sources and usages ~~",
             "}"),
             examples = c("\\examples{",
             paste0("data(", name, ")"),
             paste0("## maybe str(", name, ") ; plot(", name, ") ..."),
             "}"),
             keywords = "\\keyword{datasets}")

    if(is.na(filename)) return(Rdtxt)

    cat(unlist(Rdtxt), file = filename, sep = "\n")

    message(gettextf("Created file named %s.", sQuote(filename)),
            "\n",
            gettext("Edit the file and move it to the appropriate directory."),
            domain = NA)

    invisible(filename)
}

promptPackage <-
function(package, lib.loc = NULL, filename = NULL, name = NULL, final = FALSE)
{
    ## Most of this should not be translated -- PR#11191
    ## need to do this as packageDescription and library(help=) have
    ## different conventions
    if (is.null(lib.loc)) lib.loc <- .libPaths()

    insert1 <- function(field, new) {
    	prev <- Rdtxt[[field]]
    	Rdtxt[[field]] <<- c(prev[-length(prev)], new, prev[length(prev)])
    }
    insert2 <- function(field, new) insert1(field, paste("~~", new, "~~"))
    tabular <- function(col1, col2)
        c("\\tabular{ll}{", paste0(col1, " \\tab ", col2, "\\cr"), "}")

    if(missing(name))
	name <- paste0(package, "-package")

    if(is.null(filename))
        filename <- paste0(name, ".Rd")

    Rdtxt <-
    	    list(name = paste0("\\name{", name, "}"),
#                 version = "\\Rdversion{1.1}",
    	         aliases = paste0("\\alias{", name, "}"),
    	         docType = "\\docType{package}",
    	         title = c("\\title{", "}"),
    	         description = c("\\description{","}"),
    	         details = c("\\details{","}"),
    	         author = c("\\author{","}"),
    	         references = character(0L),

    	         keywords = c("\\keyword{ package }")
    	     )

    desc <- packageDescription(package, lib.loc)

    if (length(desc) > 1) {
    	info <- library(help = package, lib.loc = lib.loc,
                        character.only = TRUE)

    	if (!length(grep(paste0("^", package, " "), info$info[[2L]])))
    	    Rdtxt$aliases <- c(Rdtxt$aliases, paste0("\\alias{", package, "}"))

        insert1("title", desc$Title)
	insert1("description", desc$Description)
	insert1("author", c(desc$Author, "",
                            paste(identity("Maintainer:"),desc$Maintainer)))

	desc <- desc[!(names(desc) %in%
                       c("Title", "Description", "Author", "Maintainer"))]

	insert1("details", tabular(paste0(names(desc), ":"), unlist(desc)))

	if (!is.null(info$info[[2L]]))
	    insert1("details",  c("", identity("Index:"), "\\preformatted{",
	                          info$info[[2L]], "}"))
	if (!is.null(info$info[[3L]]))
	    insert1("details",
                    c("",
        identity("Further information is available in the following vignettes:"),
                      tabular(paste0("\\code{", info$info[[3L]][,1], "}"),
                              info$info[[3L]][,2])))
    }

    if (!final) {
        insert2("title", identity("package title"))
        insert2("description",
                identity("A concise (1-5 lines) description of the package"))
        insert2("details",
                strwrap(identity("An overview of how to use the package, including the most important functions")))
        insert2("author",
                identity("The author and/or maintainer of the package"))
        Rdtxt$references <-
            c("\\references{",
              paste("~~",
                    identity("Literature or other references for background information"),
                    "~~"),
              "}")
        Rdtxt$seealso <- c("\\seealso{", "}")
        insert2("seealso",
                c(identity("Optional links to other man pages, e.g."),
                  "\\code{\\link[<pkg>:<pkg>-package]{<pkg>}}"))
        Rdtxt$examples <- c("\\examples{","}")
        insert2("examples",
                identity("simple examples of the most important functions"))
        insert2("keywords",
                strwrap(identity("Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory")))
    }

    if(is.na(filename)) return(Rdtxt)

    cat(unlist(Rdtxt), file = filename, sep = "\n")

    message(gettextf("Created file named %s.", sQuote(filename)),
            "\n",
            gettext("Edit the file and move it to the appropriate directory."),
            domain = NA)

    invisible(filename)
}
#  File src/library/utils/R/question.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

`?` <-
function(e1, e2)
{
    if (missing(e2)) {
	type <- NULL
	topicExpr <- substitute(e1)
    } else {
	type <- substitute(e1)
	topicExpr <- substitute(e2)
    }
    search <- (is.call(topicExpr) && topicExpr[[1L]] == "?")
    if(search) { # ??foo is parsed as `?`(`?`(foo))
	topicExpr <- topicExpr[[2L]]
	if (is.call(te <- topicExpr	 ) && te[[1L]] == "?" &&
	    is.call(te <- topicExpr[[2L]]) && te[[1L]] == "?") {
	    cat("Contacting Delphi...")
	    flush.console()
	    Sys.sleep(2+rpois(1,2))
	    cat("the oracle is unavailable.\nWe apologize for any inconvenience.\n")
	    return(invisible())
	}
    }

    if (is.call(topicExpr) && (topicExpr[[1L]] == "::" ||
			       topicExpr[[1L]] == ":::")) {
	package <- as.character(topicExpr[[2L]])
	topicExpr <- topicExpr[[3L]]
    }
    else
	package <- NULL

    if (search) {
	if(is.null(type))
	    return(eval(substitute(help.search(TOPIC, package = PACKAGE),
				   list(TOPIC = as.character(topicExpr),
					PACKAGE = package))))
	else
	    return(eval(substitute(help.search(TOPIC, fields = FIELD,
					       package = PACKAGE),
				   list(TOPIC = as.character(topicExpr),
					FIELD = as.character(type),
					PACKAGE = package))))
    } else {
	if (is.null(type)) {
	    if (is.call(topicExpr))
		return(.helpForCall(topicExpr, parent.frame()))
	    topic <-
		if(is.name(topicExpr)) as.character(topicExpr) else e1
	    return(eval(substitute(help(TOPIC, package = PACKAGE),
				   list(TOPIC = topic,
					PACKAGE = package))))
	} else {
	    ## interpret e1 as a type, but to allow customization, do NOT
	    ## force arbitrary expressions to be single character strings
	    ## (so that methods can be defined for topicName).
	    type <-
		if(is.name(type)) as.character(type) else e1
	    topic <-
		if(is.name(topicExpr)) as.character(topicExpr)
		else {
		    if (is.call(topicExpr) && identical(type, "method"))
			return(.helpForCall(topicExpr, parent.frame(), FALSE))
		    e2
		}
            h <- .tryHelp(topicName(type, topic), package = package)
            if(is.null(h)) {
		if(is.language(topicExpr))
		    topicExpr <- deparse(topicExpr)
		stop(gettextf("no documentation of type %s and topic %s (or error in processing help)",
			      sQuote(type), sQuote(topicExpr)),
                     domain = NA)
	    }
            h
	}
    }
}

topicName <-
function(type, topic)
{
    if((length(type) == 0L) || (length(topic) == 0L))
        character(0L)
    else
        paste(paste(topic, collapse = ","), type, sep = "-")
}

.helpForCall <-
function(expr, envir, doEval = TRUE)
{
    ## There should really be a common way of formatting signatures.
    sigFormat <- function(sigNames, sigClasses) {
        paste(sprintf("%s = \"%s\"", sigNames, sigClasses),
              collapse = ", ")
    }

    f <- expr[[1L]]                     # the function specifier
    where <- topenv(envir)              # typically .GlobalEnv
    if(is.name(f))
        f <- as.character(f)
    if(!.isMethodsDispatchOn() || !methods::isGeneric(f, where = where)) {
        if(!is.character(f) || length(f) != 1L)
            stop(gettextf("the object of class %s in the function call %s could not be used as a documentation topic",
                          dQuote(class(f)), sQuote(deparse(expr))),
                 domain = NA)
        h <- .tryHelp(f)
        if(is.null(h))
            stop(gettextf("no methods for %s and no documentation for it as a function",
                          sQuote(f)),
                 domain = NA)
    }
    else {
        ## allow generic function objects or names
        if(methods::is(f, "genericFunction")) {
            fdef <- f
            f <- fdef@generic
        }
        else
            fdef <- methods::getGeneric(f, where = where)
        call <- match.call(fdef, expr)
        ## make the signature
        sigNames <- fdef@signature
        sigClasses <- rep.int("ANY", length(sigNames))
        names(sigClasses) <- sigNames
        for(arg in sigNames) {
            argExpr <- methods::elNamed(call, arg)
            if(!is.null(argExpr)) {
                simple <- (is.character(argExpr) || is.name(argExpr))
                ## TODO:  ideally, if doEval is TRUE, we would like to
                ## create the same context used by applyClosure in
                ## eval.c, but then skip the actual evaluation of the
                ## body.  If we could create this environment then
                ## passing it to selectMethod is closer to the semantics
                ## of the "real" function call than the code below.
                ## But, seems to need a change to eval.c and a flag to
                ## the evaluator.
                if(doEval || !simple) {
                    argVal <- try(eval(argExpr, envir))
                    if(methods::is(argVal, "try-error"))
                        stop(gettextf("error in trying to evaluate the expression for argument %s (%s)",
                                      sQuote(arg), deparse(argExpr)),
                             domain = NA)
                    sigClasses[[arg]] <- class(argVal)
                }
                else
                    sigClasses[[arg]] <- as.character(argExpr)
            }
        }
        method <- methods::selectMethod(f, sigClasses, optional=TRUE,
                                        fdef = fdef)
        if(methods::is(method, "MethodDefinition")) {
            sigClasses <- method@defined
            if(length(sigClasses) < length(sigNames))
                sigClasses <-
                    c(sigClasses,
                      rep.int("ANY", length(sigNames) - length(sigClasses)))
        }
        else
            warning(gettextf("no method defined for function %s and signature %s",
                             sQuote(f),
                             sQuote(sigFormat(sigNames, sigClasses))),
                    domain = NA)
        topic <- topicName("method", c(f, sigClasses))
        h <- .tryHelp(topic)
        if(is.null(h))
            stop(gettextf("no documentation for function %s and signature %s",
                          sQuote(f),
                          sQuote(sigFormat(sigNames, sigClasses))),
                 domain = NA)
    }

    h
}

.tryHelp <-
function(topic, package = NULL)
{
    ## Try finding help.
    ## Return NULL (nothing) in case we found no help pages, or an
    ## error.
    ## (Earlier versions showed what they found via print(), or gave
    ## an error.)
    h <- tryCatch(do.call("help", list(topic, package = package)),
                  error = identity)
    if(inherits(h, "error") || !length(h))
        h <- NULL
    h
}
#  File src/library/utils/R/read.DIF.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

read.DIF <- function(file, header = FALSE, dec = ".",
         row.names, col.names, as.is = !stringsAsFactors,
         na.strings = "NA", colClasses = NA,
         nrows = -1, skip = 0,
         check.names = TRUE,
         blank.lines.skip = TRUE,
         stringsAsFactors = default.stringsAsFactors(),
	 transpose = FALSE)
{
    if (.Platform$OS.type == "windows" && identical(file, "clipboard")) {
	if (!(5 %in% getClipboardFormats(numeric=TRUE)) ) stop("No DIF data on clipboard")
	lines <- readClipboard(5)
    } else {
	lines <- readLines(file)
    }
    if(length(lines) < 1L) stop("file had no lines")
    topic <- ""
    nrow <- NA
    ncol <- NA
    i <- 1L
    ## Read header info :
    while (topic != "DATA") {
	topic <- lines[i]
	vnum <- lines[i+1]
	num <- as.numeric(sub("^.*,","",vnum))
	## v <- as.numeric(sub(",.*$","",vnum))
	## value <- lines[i+2]
	i <- i + 3L
	if (topic == "VECTORS")
	    if(transpose) nrow <- num else ncol <- num
	else if (topic == "TUPLES")
	    if(transpose) ncol <- num else nrow <- num
    }
    if (is.na(nrow) || is.na(ncol)) stop("row and column counts not found")

    data <- matrix("", nrow, ncol)
    types <- matrix(NA_character_, nrow, ncol)

    row <- 0L
    while (i < length(lines)) {
	typenum <- lines[i]
	type <- as.numeric(sub(",.*$","",typenum))
	num <- as.numeric(sub("^.*,","",typenum))
	stringval <- lines[i+1]
	i <- i + 2L
	if (type == -1L) {
	    if (stringval == "BOT") {
		row <- row + 1L
                if(row > nrow)
                    stop("More rows than specified in header; maybe use 'transpose=TRUE'")
		col <- 0L
	    } else if (stringval == "EOD") break
	    else stop("Unrecognized special data value")
	} else {
	    col <- col + 1L
            if(col > ncol)
                stop("More columns than specified in header; maybe use 'transpose=TRUE'")
            if (type == 0L) {
                types[row, col] <- "numeric"
                if (stringval == "V") data[row, col] <- num
                else if (stringval == "NA") data[row, col] <- NA
                else if (stringval == "ERROR") data[row, col] <- NA
                else if (stringval == "TRUE") {
                    data[row, col] <- "TRUE"
                    types[row, col] <- "logical"
                }
                else if (stringval == "FALSE") {
                    data[row, col] <- "FALSE"
                    types[row, col] <- "logical"
                }
                else stop("Unrecognized value indicator")
            } else if (type == 1L) {
                types[row, col] <- "character"
                stringval <- sub("^\"", "", stringval)
                stringval <- sub("\"$", "", stringval)
                data[row, col] <- stringval
            }
        }
    }

    if(skip > 0L) data <- data[-(1L:skip),,drop=FALSE]

    ## determine header, no of cols.
    nlines <- nrow(data)

    if (!nlines) {
        if (missing(col.names))
            stop("no lines available in input")
        else {
            tmp <- vector("list", length(col.names))
            names(tmp) <- col.names
            class(tmp) <- "data.frame"
            return(tmp)
        }
    }
    first <- data[1L, ]
    if (first[1L] == "") first <- first[-1L]

    col1 <- if(missing(col.names)) length(first) else length(col.names)
    cols <- ncol

    ##	basic column counting and header determination;
    ##	rlabp (logical) := it looks like we have column names
    rlabp <- all(types[1L, ][-1L] == "character") && data[1L, 1L] == ""
    if(rlabp && missing(header))
	header <- TRUE
    if(!header) rlabp <- FALSE

    if (header) {
    	data <- data[-1L,,drop=FALSE] # skip over header
    	types <- types[-1L,,drop=FALSE]
        if(missing(col.names)) col.names <- first
        else if(length(first) != length(col.names))
            warning("header and 'col.names' are of different lengths")

    } else if (missing(col.names))
	col.names <- paste0("V", 1L:cols)
    if(length(col.names) + rlabp < cols)
        stop("more columns than column names")
    if(cols > 0L && length(col.names) > cols)
        stop("more column names than columns")
    if(cols == 0L) stop("rows are empty: giving up")


    if(check.names) col.names <- make.names(col.names, unique = TRUE)
    if (rlabp) col.names <- c("row.names", col.names)

    nmColClasses <- names(colClasses)
    if(length(colClasses) < cols)
        if(is.null(nmColClasses)) {
            colClasses <- rep_len(colClasses, cols)
        } else {
            tmp <- rep_len(NA_character_, cols)
            names(tmp) <- col.names
            i <- match(nmColClasses, col.names, 0L)
            if(any(i <= 0L))
                warning("not all columns named in 'colClasses' exist")
            tmp[ i[i > 0L] ] <- colClasses
            colClasses <- tmp
        }


    ##	set up as if we'll scan the file.

    colClasses[colClasses %in% c("real", "double")] <- "numeric"
    known <- colClasses %in%
                c("logical", "integer", "numeric", "complex", "character")
    keep <- !(colClasses %in% "NULL")

    if (blank.lines.skip) data <- data[apply(data, 1L, function(x) !all(x == "")),,drop=FALSE]
    if (nrows > -1 && nrows < nrow(data)) data <- data[seq_len(nrows),,drop=FALSE]
    nlines <- nrow(data)

    data[data %in% na.strings] <- NA
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    names(data) <- col.names

    ##	now we have the data;
    ##	convert to numeric or factor variables
    ##	(depending on the specified value of "as.is").
    ##	we do this here so that columns match up

    if(cols != length(data)) { # this should never happen
	warning("cols = ", cols, " != length(data) = ", length(data),
                domain = NA)
	cols <- length(data)
    }

    if(is.logical(as.is)) {
	as.is <- rep_len(as.is, cols)
    } else if(is.numeric(as.is)) {
	if(any(as.is < 1 | as.is > cols))
	    stop("invalid numeric 'as.is' expression")
	i <- rep.int(FALSE, cols)
	i[as.is] <- TRUE
	as.is <- i
    } else if(is.character(as.is)) {
        i <- match(as.is, col.names, 0L)
        if(any(i <= 0L))
            warning("not all columns named in 'as.is' exist")
        i <- i[i > 0L]
        as.is <- rep.int(FALSE, cols)
        as.is[i] <- TRUE
    } else if (length(as.is) != cols)
	stop(gettextf("'as.is' has the wrong length %d  != cols = %d",
                     length(as.is), cols), domain = NA)

    do <- keep & !known # & !as.is
    if(rlabp) do[1L] <- FALSE # don't convert "row.names"
    for (i in (1L:cols)[do]) {
        data[[i]] <-
	    if (is.na(colClasses[i])) {
	        if (any(types[,i] == "character")) {
	            if (stringsAsFactors && !as.is[i]) as.factor(data[[i]])
	            else data[[i]]
		} else
		    type.convert(data[[i]], as.is = as.is[i], dec = dec,
				 na.strings = character(0L))
	    }
        ## as na.strings have already been converted to <NA>
            else if (colClasses[i] == "factor") as.factor(data[[i]])
            else if (colClasses[i] == "Date") as.Date(data[[i]])
            else if (colClasses[i] == "POSIXct") as.POSIXct(data[[i]])
            else methods::as(data[[i]], colClasses[i])
    }

    ##	now determine row names
    compactRN <- TRUE
    if (missing(row.names)) {
	if (rlabp) {
	    row.names <- data[[1L]]
	    data <- data[-1L]
            keep <- keep[-1L]
            compactRN <- FALSE
	}
	else row.names <- .set_row_names(as.integer(nlines))
    } else if (is.null(row.names)) {
	row.names <- .set_row_names(as.integer(nlines))
    } else if (is.character(row.names)) {
        compactRN <- FALSE
	if (length(row.names) == 1L) {
	    rowvar <- (1L:cols)[match(col.names, row.names, 0L) == 1L]
	    row.names <- data[[rowvar]]
	    data <- data[-rowvar]
            keep <- keep[-rowvar]
	}
    } else if (is.numeric(row.names) && length(row.names) == 1L) {
        compactRN <- FALSE
	rlabp <- row.names
	row.names <- data[[rlabp]]
	data <- data[-rlabp]
        keep <- keep[-rlabp]
    } else stop("invalid 'row.names' specification")
    data <- data[keep]

    ## rownames<- is interpreted, so avoid it for efficiency (it will copy)
    if(is.object(row.names) || !(is.integer(row.names)) )
        row.names <- as.character(row.names)
    if(!compactRN) {
        if (length(row.names) != nlines)
            stop("invalid 'row.names' length")
        if (anyDuplicated(row.names))
            stop("duplicate 'row.names' are not allowed")
        if (any(is.na(row.names)))
            stop("missing values in 'row.names' are not allowed")
    }

    ##	this is extremely underhanded
    ##	we should use the constructor function ...
    ##	don't try this at home kids

    attr(data, "row.names") <- row.names
    data
}
#  File src/library/utils/R/read.fortran.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


read.fortran <- function(file, format, ..., as.is = TRUE, colClasses = NA)
{

    processFormat <- function(format){
        format <- toupper(format)
        template <- "^([0-9]*)([FXAI])([0-9]*)\\.?([0-9]*)"
        reps <- as.numeric(sub(template,"\\1",format))
        types <- sub(template, "\\2", format)
        lengths <- as.numeric(sub(template, "\\3", format))
        decimals <- as.numeric(sub(template, "\\4", format))

        reps[is.na(reps)] <- 1L
        lengths[is.na(lengths) & types=="X"] <- 1L

        charskip <- types=="X"
        lengths[charskip] <- reps[charskip]*lengths[charskip]
        reps[charskip] <- 1

        if (any(is.na(lengths)))
            stop("missing lengths for some fields")

        lengths <- rep(lengths,reps)
        types <- rep(types,reps)
        decimals <- rep(decimals,reps)
        types <-  match(types, c("F","D","X","A","I"))

        if (any(!is.na(decimals) & types>2L))
            stop("invalid format")
        colClasses  <-  c("numeric", "numeric", NA,
                          if(as.is) "character" else NA, "integer")[types]
        colClasses  <-  colClasses[!(types==3L)]
        decimals  <-    decimals  [!(types==3L)]
        lengths[types==3] <-  -lengths[types==3L]

        list(lengths,colClasses,decimals)
    }

    if(is.list(format)){
        ff <- lapply(format,processFormat)
        widths <- lapply(ff,"[[",1L)
        if (is.na(colClasses))
            colClasses <- do.call("c",lapply(ff,"[[",2L))
        decimals <- do.call("c",lapply(ff,"[[",3L))
    } else {
        ff <- processFormat(format)
        widths <- ff[[1L]]
        if (is.na(colClasses))
            colClasses <- ff[[2L]]
        decimals <- ff[[3L]]
    }
    rval <- read.fwf(file,widths=widths, ..., colClasses=colClasses)
    for(i in which(!is.na(decimals)))
        rval[,i] <- rval[,i]*(10^-decimals[i])
    rval
}
#  File src/library/utils/R/read.fwf.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

read.fwf <-
function(file, widths, header = FALSE, sep = "\t",
         skip = 0, row.names, col.names, n = -1, buffersize = 2000, ...)
{
    doone <- function(x) {
        x <- substring(x, first, last)
        x[!nzchar(x)] <- NA_character_
        x
    }

    if (is.list(widths)) {
        recordlength <- length(widths)
        widths <- do.call("c", widths)
    } else recordlength <- 1L

    drop <- (widths < 0L)
    widths <- abs(widths)


    buffersize <- (buffersize %/% recordlength) * recordlength

    FILENAME <- tempfile("Rfwf.")
    on.exit(unlink(FILENAME))
    FILE <- file(FILENAME,"a")
    on.exit(close(FILE),add=TRUE)

    if (is.character(file)) {
        file <- file(file, "rt")
        on.exit(close(file), add=TRUE)
    } else if (!isOpen(file)) {
        open(file, "rt")
        on.exit(close(file), add=TRUE)
    }

    if (skip) readLines(file, n=skip)
    if (header) {
        headerline <- readLines(file, n=1L)
        cat(file=FILE, headerline, "\n")
    }

    repeat({
        if (n == 0L) break
        if (n == -1L)
            thisblock <- buffersize
        else
            thisblock <- min(buffersize,n*recordlength)

        raw <- readLines(file, n = thisblock)
        nread <- length(raw)
        if (recordlength > 1L &&  nread %% recordlength) {
            raw <- raw[1L:(nread-nread %% recordlength)]
            warning(sprintf(ngettext(nread %% recordlength,
                                     "last record incomplete, %d line discarded",
                                     "last record incomplete, %d lines discarded"),
                            nread %% recordlength), domain = NA)
        }
        if (recordlength > 1L) {
            raw <- matrix(raw, nrow=recordlength)
            raw <- apply(raw, 2L, paste, collapse="")
        }

        st <- c(1L, 1L+cumsum(widths))
        first <- st[-length(st)][!drop]
        last <- cumsum(widths)[!drop]
        cat(file = FILE, sapply(raw, doone),
            sep = c(rep_len(sep, length(first)-1L), "\n"))

        if (nread < thisblock) break
        if (n > 0L) n <- n - length(raw)
    })

    close(FILE)
    FILE <- file(FILENAME,"r")
    read.table(file = FILE, header = header, sep = sep,
	       row.names = row.names, col.names = col.names, quote = "", ...)
}
#  File src/library/utils/R/readhttp.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

url.show <-
    function (url,  title = url, file = tempfile(),
              delete.file = TRUE, method, ...)
{
    if (download.file(url, destfile = file, method = method, mode = "w"))
        stop("transfer failure")
    file.show(file, delete.file = delete.file, title = title, ...)
}



defaultUserAgent <- function()
{
    Rver <- paste(R.version$major, R.version$minor, sep=".")
    Rdetails <- paste(Rver, R.version$platform, R.version$arch,
                      R.version$os)
    paste0("R (", Rdetails, ")")
}


makeUserAgent <- function(format = TRUE) {
    agent <- getOption("HTTPUserAgent")
    if (is.null(agent)) {
        return(NULL)
    }
    if (length(agent) != 1L)
        stop(gettextf("%s option must be a length one character vector or NULL",
                      sQuote("HTTPUserAgent")), domain = NA)
    if (format) paste0("User-Agent: ", agent[1L], "\r\n") else agent[1L]
}
#  File src/library/utils/R/readtable.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

count.fields <-
function(file, sep = "", quote = "\"'", skip = 0,
         blank.lines.skip = TRUE, comment.char = "#")
{
    if(is.character(file)) {
        file <- file(file)
        on.exit(close(file))
    }
    if(!inherits(file, "connection"))
        stop("'file' must be a character string or connection")
    .External(C_countfields, file, sep, quote, skip, blank.lines.skip,
              comment.char)
}


type.convert <-
function(x, na.strings = "NA", as.is = FALSE, dec = ".")
    .External2(C_typeconvert, x, na.strings, as.is, dec)


read.table <-
function(file, header = FALSE, sep = "", quote = "\"'", dec = ".",
         row.names, col.names, as.is = !stringsAsFactors,
         na.strings = "NA", colClasses = NA,
         nrows = -1, skip = 0,
         check.names = TRUE, fill = !blank.lines.skip,
         strip.white = FALSE, blank.lines.skip = TRUE,
         comment.char = "#", allowEscapes = FALSE, flush = FALSE,
         stringsAsFactors = default.stringsAsFactors(),
         fileEncoding = "", encoding = "unknown", text)
{
    if (missing(file) && !missing(text)) {
	file <- textConnection(text)
	on.exit(close(file))
    }
    if(is.character(file)) {
        file <- if(nzchar(fileEncoding))
            file(file, "rt", encoding = fileEncoding) else file(file, "rt")
        on.exit(close(file))
    }
    if(!inherits(file, "connection"))
        stop("'file' must be a character string or connection")
    if(!isOpen(file, "rt")) {
        open(file, "rt")
        on.exit(close(file))
    }

    if(skip > 0L) readLines(file, skip)
    ## read a few lines to determine header, no of cols.
    nlines <- n0lines <- if (nrows < 0L) 5 else min(5L, (header + nrows))

    lines <- .External(C_readtablehead, file, nlines, comment.char,
                       blank.lines.skip, quote, sep)
    nlines <- length(lines)
    if(!nlines) {
        if(missing(col.names)) stop("no lines available in input")
        rlabp <- FALSE
        cols <- length(col.names)
    } else {
        if(all(!nzchar(lines))) stop("empty beginning of file")
        if(nlines < n0lines && file == 0L)  { # stdin() has reached EOF
            pushBack(c(lines, lines, ""), file)
            on.exit((clearPushBack(stdin())))
        } else pushBack(c(lines, lines), file)
        first <- scan(file, what = "", sep = sep, quote = quote,
                      nlines = 1, quiet = TRUE, skip = 0,
                      strip.white = TRUE,
                      blank.lines.skip = blank.lines.skip,
                      comment.char = comment.char, allowEscapes = allowEscapes,
                      encoding = encoding)
        col1 <- if(missing(col.names)) length(first) else length(col.names)
        col <- numeric(nlines - 1L)
        if (nlines > 1L)
            for (i in seq_along(col))
                col[i] <- length(scan(file, what = "", sep = sep,
                                      quote = quote,
                                      nlines = 1, quiet = TRUE, skip = 0,
                                      strip.white = strip.white,
                                      blank.lines.skip = blank.lines.skip,
                                      comment.char = comment.char,
                                      allowEscapes = allowEscapes))
        cols <- max(col1, col)

        ##	basic column counting and header determination;
        ##	rlabp (logical) := it looks like we have column names

        rlabp <- (cols - col1) == 1L
        if(rlabp && missing(header))
            header <- TRUE
        if(!header) rlabp <- FALSE

        if (header) {
            ## skip over header
           .External(C_readtablehead, file, 1L, comment.char,
                     blank.lines.skip, quote, sep)
            if(missing(col.names)) col.names <- first
            else if(length(first) != length(col.names))
                warning("header and 'col.names' are of different lengths")

        } else if (missing(col.names))
            col.names <- paste0("V", 1L:cols)
        if(length(col.names) + rlabp < cols)
            stop("more columns than column names")
        if(fill && length(col.names) > cols)
            cols <- length(col.names)
        if(!fill && cols > 0L && length(col.names) > cols)
            stop("more column names than columns")
        if(cols == 0L) stop("first five rows are empty: giving up")
    }

    if(check.names) col.names <- make.names(col.names, unique = TRUE)
    if (rlabp) col.names <- c("row.names", col.names)

    nmColClasses <- names(colClasses)
    if(length(colClasses) < cols)
        if(is.null(nmColClasses)) {
            colClasses <- rep_len(colClasses, cols)
        } else {
            tmp <- rep_len(NA_character_, cols)
            names(tmp) <- col.names
            i <- match(nmColClasses, col.names, 0L)
            if(any(i <= 0L))
                warning("not all columns named in 'colClasses' exist")
            tmp[ i[i > 0L] ] <- colClasses
            colClasses <- tmp
        }


    ##	set up for the scan of the file.
    ##	we read unknown values as character strings and convert later.

    what <- rep.int(list(""), cols)
    names(what) <- col.names

    colClasses[colClasses %in% c("real", "double")] <- "numeric"
    known <- colClasses %in% c("logical", "integer", "numeric", "complex",
                               "character", "raw")
    what[known] <- sapply(colClasses[known], do.call, list(0))
    what[colClasses %in% "NULL"] <- list(NULL)
    keep <- !sapply(what, is.null)

    data <- scan(file = file, what = what, sep = sep, quote = quote,
                 dec = dec, nmax = nrows, skip = 0,
		 na.strings = na.strings, quiet = TRUE, fill = fill,
                 strip.white = strip.white,
                 blank.lines.skip = blank.lines.skip, multi.line = FALSE,
                 comment.char = comment.char, allowEscapes = allowEscapes,
                 flush = flush, encoding = encoding)

    nlines <- length(data[[ which.max(keep) ]])

    ##	now we have the data;
    ##	convert to numeric or factor variables
    ##	(depending on the specified value of "as.is").
    ##	we do this here so that columns match up

    if(cols != length(data)) { # this should never happen
	warning("cols = ", cols, " != length(data) = ", length(data),
                domain = NA)
	cols <- length(data)
    }

    if(is.logical(as.is)) {
	as.is <- rep_len(as.is, cols)
    } else if(is.numeric(as.is)) {
	if(any(as.is < 1 | as.is > cols))
	    stop("invalid numeric 'as.is' expression")
	i <- rep.int(FALSE, cols)
	i[as.is] <- TRUE
	as.is <- i
    } else if(is.character(as.is)) {
        i <- match(as.is, col.names, 0L)
        if(any(i <= 0L))
            warning("not all columns named in 'as.is' exist")
        i <- i[i > 0L]
        as.is <- rep.int(FALSE, cols)
        as.is[i] <- TRUE
    } else if (length(as.is) != cols)
	stop(gettextf("'as.is' has the wrong length %d  != cols = %d",
                     length(as.is), cols), domain = NA)

    do <- keep & !known # & !as.is
    if(rlabp) do[1L] <- FALSE # don't convert "row.names"
    for (i in (1L:cols)[do]) {
        data[[i]] <-
            if (is.na(colClasses[i]))
                type.convert(data[[i]], as.is = as.is[i], dec = dec,
                             na.strings = character(0L))
        ## as na.strings have already been converted to <NA>
            else if (colClasses[i] == "factor") as.factor(data[[i]])
            else if (colClasses[i] == "Date") as.Date(data[[i]])
            else if (colClasses[i] == "POSIXct") as.POSIXct(data[[i]])
            else methods::as(data[[i]], colClasses[i])
    }

    ##	now determine row names
    compactRN <- TRUE
    if (missing(row.names)) {
	if (rlabp) {
	    row.names <- data[[1L]]
	    data <- data[-1L]
            keep <- keep[-1L]
            compactRN <- FALSE
	}
	else row.names <- .set_row_names(as.integer(nlines))
    } else if (is.null(row.names)) {
	row.names <- .set_row_names(as.integer(nlines))
    } else if (is.character(row.names)) {
        compactRN <- FALSE
	if (length(row.names) == 1L) {
	    rowvar <- (1L:cols)[match(col.names, row.names, 0L) == 1L]
	    row.names <- data[[rowvar]]
	    data <- data[-rowvar]
            keep <- keep[-rowvar]
	}
    } else if (is.numeric(row.names) && length(row.names) == 1L) {
        compactRN <- FALSE
	rlabp <- row.names
	row.names <- data[[rlabp]]
	data <- data[-rlabp]
        keep <- keep[-rlabp]
    } else stop("invalid 'row.names' specification")
    data <- data[keep]

    ## rownames<- is interpreted, so avoid it for efficiency (it will copy)
    if(is.object(row.names) || !(is.integer(row.names)) )
        row.names <- as.character(row.names)
    if(!compactRN) {
        if (length(row.names) != nlines)
            stop("invalid 'row.names' length")
        if (anyDuplicated(row.names))
            stop("duplicate 'row.names' are not allowed")
        if (any(is.na(row.names)))
            stop("missing values in 'row.names' are not allowed")
    }

    ##	this is extremely underhanded
    ##	we should use the constructor function ...
    ##	don't try this at home kids

    class(data) <- "data.frame"
    attr(data, "row.names") <- row.names
    data
}

read.csv <-
function (file, header = TRUE, sep = ",", quote = "\"", dec = ".",
          fill = TRUE, comment.char = "", ...)
    read.table(file = file, header = header, sep = sep,
               quote = quote, dec = dec, fill = fill,
               comment.char = comment.char,  ...)

read.csv2 <-
function (file, header = TRUE, sep = ";", quote = "\"", dec = ",",
          fill = TRUE, comment.char = "", ...)
    read.table(file = file, header = header, sep = sep,
               quote = quote, dec = dec, fill = fill,
               comment.char = comment.char, ...)

read.delim <-
function (file, header = TRUE, sep = "\t", quote = "\"", dec = ".",
          fill = TRUE, comment.char = "", ...)
    read.table(file = file, header = header, sep = sep,
               quote = quote, dec = dec, fill = fill,
               comment.char = comment.char, ...)

read.delim2 <-
function (file, header = TRUE, sep = "\t", quote = "\"", dec = ",",
          fill = TRUE, comment.char = "", ...)
    read.table(file = file, header = header, sep = sep,
               quote = quote, dec = dec, fill = fill,
               comment.char = comment.char, ...)

#  File src/library/utils/R/relist.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

# relist.R -- an inverse operator to unlist
# written by Andrew Clausen <clausen@econ.upenn.edu> in 2007
# with helpful suggestions from
#	Martin Maechler  ETHZ CH
#	Gabor Grothendieck at Gmail dotcom
#	Seth Falcon, sfalcon (near) fhcrc (comma) org
#
# Some functions need many parameters, which are most easily represented in
# complex structures.  Unfortunately, many mathematical functions in R,
# including optim, nlm, and grad can only operate on functions whose domain is
# a vector.  R has a function called "unlist" to convert complex objects into a
# vector representation.  This file provides an inverse operation called
# "relist" to convert vectors back to the convenient structural representation.
# Together, these functions allow structured functions to have simple
# mathematical interfaces.
#
# For example, a likelihood function for a multivariate normal model needs a
# variance-covariance matrix and a mean vector.	 It would be most convenient to
# represent it as a list containing a vector and a matrix.  A typical parameter
# might look like
#
#	list(mean=c(0, 1), vcov=cbind(c(1, 1), c(1, 0)))
#
# However, optim can't operate on functions that take lists as input; it
# only likes vectors.  The solution is conversion:
#
## 	initial.param <- list(mean=c(0, 1), vcov=cbind(c(1, 1), c(1, 0)))
## 	initial.param <- as.relistable(initial.param)
## #
## 	ll <- function(param.vector)
## 	{
## 		param <- relist(initial.param)
## 		-sum(dnorm(x, mean=param$mean, vcov=param$vcov, log=TRUE))
## 		# note: dnorm doesn't do vcov... but I hope you get the point
## 	}

## 	optim(unlist(initial.param), ll)
#
# "relist" takes two parameters: skeleton and flesh.  Skeleton is a sample
# object that has the right "shape" but the wrong content.  "flesh" is a vector
# with the right content but the wrong shape.  Invoking
#
#	relist(flesh, skeleton)
#
# will put the content of flesh on the skeleton.  You don't need to specify
# skeleton explicitly if the skeleton is stored as an attribute inside flesh.
# In particular, flesh was created from some object obj with
#
#	unlist(as.relistable(obj))
#
# then the skeleton attribute is automatically set.
#
# As long as "skeleton" has the right shape, it should be a precise inverse
# of unlist.  These equalities hold:
#
#	relist(unlist(x), skeleton) == x
#	unlist(relist(y, skeleton)) == y
#
#	x <- as.relistable(x)
#	relist(unlist(x)) == x

is.relistable <- function(x) inherits(x, "relistable")

as.relistable <- function(x)
{
    if (!inherits(x, "relistable"))
	class(x) <- c("relistable", class(x))
    x
}

## NB: unlist() is generic *internally* (i.e. not visible from 'unlist')
unlist.relistable <- function(x, recursive=TRUE, use.names=TRUE)
{
    if (!recursive)
	warning("relist() requires recursively unlisted objects.")
    skeleton <- x
### MM: FIXME?  I think this is just  NextMethod()
    ## remove 'relistable'
    class(x) <- setdiff(class(x), "relistable")
    result <- unlist(x, recursive, use.names)
    attr(result, "skeleton") <- skeleton
    result
}

relist <- function(flesh, skeleton=attr(flesh, "skeleton"))
{
    if (is.null(skeleton)) {
	stop("The 'flesh' argument does not contain a skeleton attribute.\n",
	     "Either ensure you unlist a relistable object, or specify the skeleton separately.")
    }
    UseMethod("relist", skeleton)
}

## was 'relist.numeric' in Andrew's code
relist.default <- function(flesh, skeleton=attr(flesh, "skeleton"))
{
    result <- flesh
    names(result) <- names(skeleton)
    result
}

relist.list <- function(flesh, skeleton=attr(flesh, "skeleton"))
{
    ind <- 1L
    result <- skeleton
    for (i in seq_along(skeleton)) {
	size <- length(unlist(result[[i]]))
	result[[i]] <-
	    relist(flesh[ind:(ind + size - 1L)], result[[i]])
	ind <- ind + size
    }
    result
}


relist.matrix <- function(flesh, skeleton=attr(flesh, "skeleton"))
{
    if (is.numeric(skeleton[1,1]))
	return(matrix(flesh, nrow=nrow(skeleton),
		      dimnames=dimnames(skeleton)))
    n <- nrow(skeleton)
    m <- ncol(skeleton)
    result <- skeleton
    ind <- 1L
    for (j in 1L:m)
	for (i in 1L:n) {
	    size <- length(unlist(skeleton[[i, j]]))
	    result[[i, j]] <- relist(flesh[ind:(ind + size - 1)],
				     skeleton[[i, j]])
	    ind <- ind + size
	}
    result
}

relist.factor <- function(flesh, skeleton=attr(flesh, "skeleton"))
{
    as.factor(levels(skeleton)[flesh])
}
#  File src/library/utils/R/roman.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

as.roman <-
function(x)
{
    if(is.numeric(x))
        x <- as.integer(x)
    else if(is.character(x)) {
        ## Let's be nice: either strings that are *all* arabics, or
        ## (hopefully, for the time being) all romans.
        x <- if(all(grepl("^[[:digit:]]+$", x)))
            as.integer(x)
        else
            .roman2numeric(x)
    }
    else
        stop("cannot coerce 'x' to roman")
    x[(x <= 0L | x >= 3900L)] <- NA
    class(x) <- "roman"
    x
}

as.character.roman <-
function(x, ...)
    .numeric2roman(x)

format.roman <-
function(x, ...)
    format(as.character(x))

print.roman <-
function(x, ...)
{
    print(noquote(as.character(x)), ...)
    invisible(x)
}

`[.roman` <-
function(x, i)
{
    cl <- oldClass(x)
    y <- NextMethod("[")
    oldClass(y) <- cl
    y
}

.numeric2roman <-
function(x) {
    romans <- c("M", "CM", "D", "CD", "C", "XC", "L", "XL", "X", "IX",
                "V", "IV", "I")
    numbers <- c(1000L, 900L, 500L, 400L, 100L, 90L, 50L, 40L, 10L, 9L,
                 5L, 4L, 1L)
    n2r <- function(z) {
        y <- character()
        for(i in seq_along(romans)) {
            d <- numbers[i]
            while(z >= d) {
                z <- z - d
                y <- c(y, romans[i])
            }
        }
        paste(y, collapse = "")
    }

    out <- character(length(x))
    x <- as.integer(x)
    ind <- is.na(x) | (x <= 0L) | (x >= 3900L)
    out[ind] <- NA
    if(any(!ind))
        out[!ind] <- sapply(x[!ind], n2r)
    out
}

.roman2numeric <-
function(x)
{
    ## <FIXME>
    ## What if this fails?
    ## Should say something like "Not a valid roman number ..."
    ## </FIXME>
    romans <- c("M", "CM", "D", "CD", "C", "XC", "L", "XL", "X", "IX",
                "V", "IV", "I")
    numbers <- c(1000L, 900L, 500L, 400L, 100L, 90L, 50L, 40L, 10L, 9L,
                 5L, 4L, 1L)
    out <- integer(length(x))
    ind <- is.na(x)
    out[ind] <- NA
    if(any(!ind)) {
        y <- toupper(x[!ind])
        y <- gsub("CM", "DCCCC", y)
        y <- gsub("CD", "CCCC", y)
        y <- gsub("XC", "LXXXX", y)
        y <- gsub("XL", "XXXX", y)
        y <- gsub("IX", "VIIII", y)
        y <- gsub("IV", "IIII", y)
        ok <- grepl("^M{,3}D?C{,4}L?X{,4}V?I{,4}$", y)
        if(any(!ok)) {
            warning(sprintf(ngettext(sum(!ok),
                                     "invalid roman numeral: %s",
                                     "invalid roman numerals: %s"),
                            paste(x[!ind][!ok], collapse = " ")),
                    domain = NA)
            out[!ind][!ok] <- NA
        }
        if(any(ok))
            out[!ind][ok] <-
                sapply(strsplit(y[ok], ""),
                       function(z)
                       as.integer(sum(numbers[match(z, romans)])))
    }
    out
}
#  File src/library/tools/R/rtags.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/



### These utilities are intended to read R code files and produce tags
### in Emacs' etags format.  Doing so in R allows us to use R's
### parser.  Support for vi-style tags could be useful, but it needs
### the tags file needs to be sorted, making file-by-file processing
### difficult. It may be easier to write a script to convert an etags
### format file (see http://http://en.wikipedia.org/wiki/Ctags).



### * shorten.to.string

## The etags format requires the initial part of a matching line to be
## recorded in the TAGS file, with an optional entry for the `token'
## that is to be matched.  Exact matches to the token are preferred,
## but it seems that subsequent non-exact matches look in the initial
## string as well as other tokens with equal priority.  Such matches
## seem pointless, so an attempt is made to shorten the matching line.
## It is not clear to me whether there are restrictions on what this
## part could be, but a completely blank string doesn't seem to work.
## For now, I'm just keeping the first letter. May change if problems
## arise.

shorten.to.string <-
    function(line, token)
{
    if (FALSE) {
        ans <- regexpr(strsplit("token", ",", fixed = TRUE)[[1L]][1L],
                       line, fixed = TRUE)
        if (ans == -1L) line
        else substr(line, 1L, ans + attr(ans, "match.length") - 1L)
    }
    else {
        ## can we just put essentially nothing? Seems to work
        substr(line, 1L, 1L)
    }
}


### * write.etags

## this function is responsible for formatting the output for a single
## file given the relevant information.  The format was inferred from
## the "Ctags" wikipedia entry and by studying etags output.

write.etags <-
    function(src,
             tokens, startlines, lines, nchars,
             ...,
             shorten.lines = c("token", "simple", "none"))
{
    ## extra 1 for newline
    shorten.lines <- match.arg(shorten.lines)
    offsets <- (cumsum(nchars + 1L) - (nchars + 1L))[startlines]
    lines <-
        switch(shorten.lines,
               none = lines,
               simple = sapply(strsplit(lines, "function", fixed = TRUE), "[", 1),
               token = mapply(shorten.to.string, lines, tokens))
    tag.lines <-
        paste(sprintf("%s\x7f%s\x01%d,%d",
                      lines, tokens, startlines,
                      as.integer(offsets)),
              collapse = "\n")
    ## simpler format: tag.lines <- paste(sprintf("%s\x7f%d,%d", lines, startlines, as.integer(offsets)), collapse = "\n")
    tagsize <- nchar(tag.lines, type = "bytes") + 1L
    cat("\x0c\n", src, ",", tagsize, "\n", tag.lines, "\n", sep = "", ...)
}


### * expr2token

## this computes the tag name from an expression.  Currently, this
## returns the second thing in the expression; so
##
##   foo <- function(x) ...      ==> `<-`,       foo, ...
##   setMethod("foo", "bar" ...  ==> setMethod,  foo, ...
##   setGeneric("foo", "bar" ... ==> setGeneric, foo, ...
##
## which covers the typical uses.  We match against a list to restrict
## types of expressions that are tagged.  To reject things like
##
##   x[i] <- 10
##
## the second component is required to have length 1.  One limitation
## is that things like
##
##   if (require(pkg)) foo <- ... else foo <- ...
##
## will not be handled.

expr2token <-
    function(x,
             ok = c("<-", "=", "<<-", "assign",
                    "setGeneric", "setGroupGeneric", "setMethod",
                    "setClass", "setClassUnion"),
             extended = TRUE)
{
    id <- ""
    value <-
        if ((length(x) > 1L) &&
            (length(token <- as.character(x[[2L]])) == 1L) &&
            (length(id <- as.character(x[[1L]])) == 1L) &&
            (id %in% ok)) token
        else
            character(0L)
    if (extended && identical(id, "setMethod"))
    {
        ## try to add the signature, comma separated
        sig <- tryCatch(eval(x[[3L]]), error = identity)
        if (!inherits(sig, "error") && is.character(sig))
            value <- paste(c(value, sig), collapse=",")
    }
    value
}


### * rtags.file

## Handles a single file

rtags.file <-
    function(src, ofile = "", append = FALSE,
             write.fun = write.etags) ## getOption("writeTags")
{

    ## FIXME: do we need to worry about encoding etc.?
    elist <- parse(src, srcfile = srcfile(src))
    if (length(elist) == 0) return(invisible())
    lines <- readLines(src)
    tokens <- lapply(elist, expr2token)
    startlines <- sapply(attr(elist, "srcref"), "[", 1L)
    if (length(tokens) != length(startlines))
        stop("length mismatch: bug in code!", domain = NA)
    keep <- sapply(tokens, length) == 1L
    if (!any(keep)) return(invisible())
    tokens <- unlist(tokens[keep])
    startlines <- startlines[keep]
    write.fun(src = src,
              tokens = tokens,
              startlines = startlines,
              lines = lines[startlines],
              nchars = nchar(lines, type = "bytes"),
              file = ofile, append = append)
}

### * rtags

## Public interface.  Tags files under a specified directory, using
## regular expressions to filter out inappropriate files.


rtags <-
    function(path = ".", pattern = "\\.[RrSs]$",
             recursive = FALSE,
             src = list.files(path = path,
                              pattern = pattern,
                              full.names = TRUE,
                              recursive = recursive),
             keep.re = NULL,
             ofile = "", append = FALSE,
             verbose = getOption("verbose"))
{
    if (ofile != "" && !append) {
        if (!file.create(ofile, showWarnings = FALSE)) 
            stop(gettextf("Could not create file %s, aborting", ofile),
                 domain = NA)
    }
    if (!missing(keep.re))
        src <- grep(keep.re, src, value = TRUE)
    for (s in src)
    {
        if (verbose) message(gettextf("Processing file %s", s), domain = NA)
        tryCatch(
                 rtags.file(s, ofile = ofile, append = TRUE),
                 error = function(e) NULL)
    }
    invisible()
}




### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***


#  File src/library/utils/R/sessionInfo.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

sessionInfo <- function(package=NULL)
{
    z <- list()
    z$R.version <- R.Version()
    z$platform <- z$R.version$platform
    if(nzchar(.Platform$r_arch))
        z$platform <- paste(z$platform, .Platform$r_arch, sep = "/")
    z$platform <- paste0(z$platform, " (", 8*.Machine$sizeof.pointer, "-bit)")
    z$locale <- Sys.getlocale()

    if(is.null(package)){
        package <- grep("^package:", search(), value=TRUE)
        # weed out environments which are not really packages
        keep <- sapply(package, function(x) x == "package:base" || !is.null(attr(as.environment(x), "path")))
        package <- sub("^package:", "", package[keep])
    }

    ## no need to re-encode given what we extract.
    pkgDesc <- lapply(package, packageDescription, encoding = NA)
    if(length(package) == 0) stop("no valid packages were specified")
    basePkgs <- sapply(pkgDesc,
                       function(x) !is.null(x$Priority) && x$Priority=="base")
    ## Hmm, see tools:::.get_standard_package_names()$base
    z$basePkgs <- package[basePkgs]
    if(any(!basePkgs)){
        z$otherPkgs <- pkgDesc[!basePkgs]
        names(z$otherPkgs) <- package[!basePkgs]
    }
    loadedOnly <- loadedNamespaces()
    loadedOnly <- loadedOnly[!(loadedOnly %in% package)]
    if (length(loadedOnly)) {
        names(loadedOnly) <- loadedOnly
        pkgDesc <- c(pkgDesc, lapply(loadedOnly, packageDescription))
        z$loadedOnly <- pkgDesc[loadedOnly]
    }
    class(z) <- "sessionInfo"
    z
}

print.sessionInfo <- function(x, locale=TRUE, ...)
{
    mkLabel <- function(L, n) {
        vers <- sapply(L[[n]], function(x) x[["Version"]])
        pkg <-  sapply(L[[n]], function(x) x[["Package"]])
        paste(pkg, vers, sep = "_")
    }

    cat(x$R.version$version.string, "\n", sep = "")
    cat("Platform: ", x$platform, "\n\n", sep = "")
    if(locale){
        cat("locale:\n")
	print(strsplit(x$locale, ";", fixed=TRUE)[[1]], quote=FALSE, ...)
        cat("\n")
    }
    cat("attached base packages:\n")
    print(x$basePkgs, quote=FALSE, ...)
    if(!is.null(x$otherPkgs)){
        cat("\nother attached packages:\n")
	print(mkLabel(x, "otherPkgs"), quote=FALSE, ...)
    }
    if(!is.null(x$loadedOnly)){
        cat("\nloaded via a namespace (and not attached):\n")
	print(mkLabel(x, "loadedOnly"), quote=FALSE, ...)
    }
    invisible(x)
}

toLatex.sessionInfo <- function(object, locale=TRUE, ...)
{
    opkgver <- sapply(object$otherPkgs, function(x) x$Version)
    nspkgver <- sapply(object$loadedOnly, function(x) x$Version)
    z <- c("\\begin{itemize}\\raggedright",
           paste0("  \\item ", object$R.version$version.string,
                  ", \\verb|", object$R.version$platform, "|"))

    if(locale){
        z <- c(z,
               paste0("  \\item Locale: \\verb|",
                      gsub(";","|, \\\\verb|", object$locale) , "|"))
    }

    z <- c(z, strwrap(paste("\\item Base packages: ",
                         paste(sort(object$basePkgs), collapse = ", ")),
                      indent = 2, exdent = 4))

    if(length(opkgver)){
        opkgver <- opkgver[sort(names(opkgver))]
        z <- c(z,
               strwrap(paste("  \\item Other packages: ",
                             paste(names(opkgver), opkgver, sep = "~",
                                   collapse = ", ")),
                       indent = 2, exdent = 4))
    }
    if(length(nspkgver)){
        nspkgver <- nspkgver[sort(names(nspkgver))]
        z <- c(z,
               strwrap(paste("  \\item Loaded via a namespace (and not attached): ",
                             paste(names(nspkgver), nspkgver, sep = "~",
                                   collapse = ", ")),
                       indent = 2, exdent = 4))
    }
    z <- c(z, "\\end{itemize}")
    class(z) <- "Latex"
    z
}
#  File src/library/utils/R/sock.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

print.socket <- function(x, ...)
{
    if(length(as.integer(x$socket)) != 1L) stop("invalid 'socket' argument")
    cat("Socket connection #", x$socket, "to", x$host, "on port", x$port, "\n")
    invisible(x)
}

make.socket <- function(host = "localhost", port, fail = TRUE, server = FALSE)
{
    if(length(port <- as.integer(port)) != 1L)
	stop("'port' must be integer of length 1")
    if(length(host <- as.character(host)) != 1L)
	stop("'host' must be character of length 1")
    if (!server){
	socket <- .Call(C_sockconnect, port, host)
    } else {
	if (host != "localhost") stop("can only receive calls on local machine")
	tmp <- .Call(C_sockopen, port)
        socket <- .Call(C_socklisten, tmp)
        host <- attr(socket, "host")
	.Call(C_sockclose, tmp)
    }
    if (socket <= 0) {
	if (fail) stop("socket not established")
        else warning("socket not established")
    }
    rval <- list(socket = socket, host, port = port)
    class(rval) <- "socket"
    rval
}

close.socket <- function(socket, ...)
    .Call(C_sockclose, socket$socket)

read.socket <- function(socket, maxlen = 256L, loop = FALSE)
{
    repeat {
	rval <- .Call(C_sockread, socket$socket, maxlen)
	if (nzchar(rval) || !loop) break
    }
    rval
}

write.socket <- function(socket, string)
    invisible(.Call(C_sockwrite, socket$socket, string))

#  File src/library/utils/R/sourceutils.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

removeSource <- function(fn) {
    stopifnot(is.function(fn))
    if (is.primitive(fn)) return(fn)
    attr(fn, "source") <- NULL
    attr(fn, "srcref") <- NULL
    attr(body(fn), "wholeSrcref") <- NULL
    attr(body(fn), "srcfile") <- NULL

    recurse <- function(part) {
        attr(part, "srcref") <- NULL
        if (is.language(part) && is.recursive(part)) {
            for (i in seq_along(part))
            	part[[i]] <- recurse(part[[i]])
        }
        part
    }
    body(fn) <- recurse(body(fn))
    fn
}

getSrcFilename <- function(x, full.names=FALSE, unique=TRUE) {
    srcref <- getSrcref(x)
    if (is.list(srcref))
    	result <- sapply(srcref, getSrcFilename, full.names, unique)
    else {
    	srcfile <- attr(srcref, "srcfile")
    	if (is.null(srcfile)) result <- character()
    	else result <- srcfile$filename
    }
    result <- if (full.names) result
              else basename(result)
    if (unique) unique(result)
    else result
}

getSrcDirectory <- function(x, unique=TRUE) {
    result <- dirname(getSrcFilename(x, full.names=TRUE, unique=unique))
    if (unique) unique(result)
    else result
}

getSrcref <- function(x) {
    if (inherits(x, "srcref")) return(x)
    if (!is.null(srcref <- attr(x, "srcref"))) return(srcref)
    if (is.function(x)) return(getSrcref(body(x)))
    NULL
}

getSrcLocation <- function(x, which=c("line", "column", "byte", "parse"), first=TRUE) {
    srcref <- getSrcref(x)
    if (is.null(srcref)) return(NULL)
    if (is.list(srcref)) sapply(srcref, getSrcLocation, which, first)
    else {
        if (length(srcref) == 6L) srcref <- c(srcref, srcref[c(1L,3L)])
    	which <- match.arg(which)
    	if (first) index <- c(line=1L, column=5L, byte=2L, parse=7L)[which]
    	else       index <- c(line=3L, column=6L, byte=4L, parse=8L)[which]
    	srcref[index]
    }
 }

getSrcfile <- function(x) {
    result <- attr(x, "srcfile")
    if (!is.null(result)) return(result)

    srcref <- attr(x, "wholeSrcref")
    if (is.null(srcref)) {
	srcref <- getSrcref(x)
    	if (is.list(srcref) && length(srcref))
    	    srcref <- srcref[[length(srcref)]]
    }
    attr(srcref, "srcfile")
}

substr_with_tabs <- function(x, start, stop, tabsize = 8) {
    widths <- rep(1, nchar(x))
    tabs <- which(strsplit(x,"")[[1]] == "\t")
    for (i in tabs) {
	cols <- cumsum(widths)
	widths[i] <- tabsize - (cols[i] - 1) %% tabsize
    }
    cols <- cumsum(widths)
    start <- which(cols >= start)
    if (!length(start))
    	return("")
    start <- start[1]
    stop <- which(cols <= stop)
    if (length(stop)) {
    	stop <- stop[length(stop)]
    	substr(x, start, stop)
    } else
    	""
}

getParseData <- function(x, includeText = NA) {
    if (inherits(x, "srcfile")) 
	srcfile <- x
    else 
	srcfile <- getSrcfile(x)

    if (is.null(srcfile))
    	return(NULL)
    else
    	data <- srcfile$parseData
    if (!is.null(data)) {
        tokens <- attr(data, "tokens")
        data <- t(unclass(data))
        colnames(data) <- c( "line1", "col1",
		 	     "line2", "col2",
		 	     "terminal", "token.num", "id", "parent" )
    	data <- data.frame(data[, -c(5,6), drop = FALSE], token = tokens,
    	                   terminal = as.logical(data[,"terminal"]),
    	                   text = attr(data, "text"),
    			   stringsAsFactors = FALSE)
    	o <- order(data[,1], data[,2], -data[,3], -data[,4])
    	data <- data[o,]
    	rownames(data) <- data$id
    	attr(data, "srcfile") <- srcfile
    	if (isTRUE(includeText)) gettext <- which(!nzchar(data$text))
        else if (is.na(includeText)) gettext <- which(!nzchar(data$text) & data$terminal)
        else {
            gettext <- integer(0)
            data$text <- NULL
        }

        if (length(gettext))
	    data$text[gettext] <- getParseText(data, data$id[gettext])
    }
    data	
}

getParseText <- function(parseData, id) {
    srcfile <- attr(parseData, "srcfile")
    d <- parseData[as.character(id),]
    text <- d$text
    if (is.null(text)) {
    	text <- character(nrow(text))
    	blank <- seq_along(text)
    } else
    	blank <- which(!nzchar(text))
    for (i in blank) {
	lines <- getSrcLines(srcfile, d$line1[i], d$line2[i])
        n <- length(lines)
        lines[n] <- substr_with_tabs(lines[n], 1, d$col2[i])
        lines[1] <- substr_with_tabs(lines[1], d$col1[i], Inf)
        text[i] <- paste(lines, collapse="\n")
    }
    text
}
#  File src/library/utils/R/str.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

####------ str : show STRucture of an R object
str <- function(object, ...) UseMethod("str")

## FIXME: convert to use sQuote
str.data.frame <- function(object, ...)
{
    ## Method to 'str' for  'data.frame' objects
    if(! is.data.frame(object)) {
	warning("str.data.frame() called with non-data.frame -- coercing to one.")
	object <- data.frame(object)
    }

    ## Show further classes // Assume that they do NOT have an own Method --
    ## not quite perfect ! (.Class = 'remaining classes', starting with current)
    cl <- oldClass(object); cl <- cl[cl != "data.frame"]  #- not THIS class
    if(0 < length(cl)) cat("Classes", paste(sQuote(cl), collapse=", "), "and ")

    cat("'data.frame':	", nrow(object), " obs. of  ",
	(p <- length(object)), " variable", if(p != 1)"s", if(p > 0)":",
	"\n", sep = "")

    ## calling next method, usually  str.default:
    if(length(l <- list(...)) && any("give.length" == names(l)))
	invisible(NextMethod("str", ...))
    else invisible(NextMethod("str", give.length=FALSE,...))
}

str.Date <- str.POSIXt <- function(object, ...) {
    cl <- oldClass(object)
    ## be careful to be fast for large object:
    n <- length(object) # FIXME, could be NA
    if(n == 0L) return(str.default(object))
    if(n > 1000L) object <- object[seq_len(1000L)]

    give.length <- TRUE ## default
    ## use 'give.length' when specified, else default = give.head
    if(length(larg <- list(...))) {
	nl <- names(larg)
	iGiveHead <- which(nl == "give.head")
	if (any(Bgl <- nl == "give.length"))
	    give.length <- larg[[which(Bgl)]]
	else if(length(iGiveHead))
	    give.length <- larg[[iGiveHead]]
	if(length(iGiveHead)) # eliminate it from arg.list
	    larg <- larg[ - iGiveHead ]
	if(is.numeric(larg[["nest.lev"]]) &&
	   is.numeric(v.len <- larg[["vec.len"]])) # typical call from data.frame
	    ## reduce length for typical call:
	    larg[["vec.len"]] <-
		min(larg[["vec.len"]],
		    (larg[["width"]]- nchar(larg[["indent.str"]]) -31)%/% 19)
    }

    le.str <- if(give.length) paste0("[1:",as.character(n),"]")
    cat(" ", cl[1L], le.str,", format: ", sep = "")
    ## do.call(str, c(list(format(object), give.head = FALSE), larg))
    ## ensuring 'object' is *not* copied:
    str.f.obj <- function(...) str(format(object), ...)
    do.call(str.f.obj, c(list(give.head = FALSE), larg))
}

strOptions <- function(strict.width = "no", digits.d = 3, vec.len = 4,
		       formatNum = function(x, ...)
		       format(x, trim=TRUE, drop0trailing=TRUE, ...))
    list(strict.width = strict.width, digits.d = digits.d, vec.len = vec.len,
	 formatNum = match.fun(formatNum))

str.default <-
    function(object, max.level = NA, vec.len = strO$vec.len,
	     digits.d = strO$digits.d,
	     nchar.max = 128, give.attr = TRUE,
	     give.head = TRUE, give.length = give.head,
	     width = getOption("width"), nest.lev = 0,
	     indent.str= paste(rep.int(" ", max(0,nest.lev+1)), collapse= ".."),
	     comp.str="$ ", no.list = FALSE, envir = baseenv(),
	     strict.width = strO$strict.width,
	     formatNum = strO$formatNum, list.len = 99,
	     ...)
{
    ## Purpose: Display STRucture of any R - object (in a compact form).
    ## --- see HELP file --
    ## ------------------------------------------------------------------------
    ## Author: Martin Maechler <maechler@stat.math.ethz.ch>	1990--1997

    ## Get defaults for these
    oDefs <- c("vec.len", "digits.d", "strict.width", "formatNum")
    ## from
    strO <- getOption("str")
    if (!is.list(strO)) {
	warning('invalid options("str") -- using defaults instead')
	strO <- strOptions()
    }
    else {
        if (!all(names(strO) %in% oDefs))
            warning(gettextf("invalid components in options(\"str\"): %s",
                             paste(setdiff(names(strO), oDefs), collapse = ", ")),
                    domain = NA)
        strO <- modifyList(strOptions(), strO)
    }
    strict.width <- match.arg(strict.width, choices = c("no", "cut", "wrap"))
    if(strict.width != "no") {
	## using eval() would be cleaner, but fails inside capture.output():
	ss <- capture.output(str.default(object, max.level = max.level,
				 vec.len = vec.len, digits.d = digits.d,
				 nchar.max = nchar.max,
				 give.attr= give.attr, give.head= give.head,
				 give.length= give.length,
				 width = width, nest.lev = nest.lev,
				 indent.str = indent.str, comp.str= comp.str,
				 no.list= no.list || is.data.frame(object),
				 envir = envir, strict.width = "no",
				 formatNum = formatNum, list.len = list.len,
					 ...) )
	if(strict.width == "wrap") {
	    nind <- nchar(indent.str) + 2
	    ss <- strwrap(ss, width = width, exdent = nind)
					# wraps at white space (only)
	}
	if(any(iLong <- nchar(ss) > width)) { ## cut hard
	    sL <- ss[iLong]
	    k <- as.integer(width-2)
	    if(any(i <- grepl("\"", substr(sL, k +1L, nchar(sL))))) {
		## care *not* to cut off the closing   "  at end of
		## string that's already truncated {-> maybe_truncate()} :
		ss[iLong[ i]] <- paste0(substr(sL[ i], 1, k-1L), "\"..")
		ss[iLong[!i]] <- paste0(substr(sL[!i], 1, k), "..")
	    } else {
		ss[iLong] <- paste0(substr(sL, 1, k),"..")
	    }
	}
	cat(ss, sep = "\n")
	return(invisible())
    }

    oo <- options(digits = digits.d); on.exit(options(oo))
    le <- length(object)
    if(is.na(le)) {
        warning("'str.default': 'le' is NA, so taken as 0", immediate. = TRUE)
        le <- 0
        vec.len <- 0
    }

    maybe_truncate <- function(x, e.x = x, Sep = "\"", ch = "| __truncated__")
    {
	trimmed <- strtrim(e.x, nchar.max)
	ii <- trimmed != e.x
	ii[is.na(ii)] <- FALSE
	if(any(ii)) x[ii] <- paste0(trimmed[ii], Sep, ch)
	x
    }
    pClass <- function(cls)
	paste0("Class", if(length(cls) > 1) "es",
	       " '", paste(cls, collapse = "', '"), "' ")
    `%w/o%` <- function(x,y) x[is.na(match(x,y))]

    nfS <- names(fStr <- formals())# names of all formal args to str.default()
    ##' Purpose: using short strSub() calls instead of long str() ones
    ##' @title Call str() on sub-parts, with mostly the *same* arguments
    ##' @param obj R object; always a "part" of the main 'object'
    ##' @param ... further arguments to str(), [often str.default()]
    strSub <- function(obj, ...) {
	## 'give.length', ...etc are *not* automatically passed down:
	nf <- nfS %w/o% c("object", "give.length", "comp.str", "no.list",
			  ## drop fn.name & "obj" :
			  names(match.call())[-(1:2)], "...")
	aList <- as.list(fStr)[nf]
	aList[] <- lapply(nf, function(n) eval(as.name(n)))
	## do.call(str, c(list(object=obj), aList, list(...)), quote=TRUE)
	## ensuring 'obj' is *not* copied:
	strObj <- function(...) str(obj, ...)
	do.call(strObj, c(aList, list(...)), quote = TRUE)
    }

    ## le.str: not used for arrays:
    le.str <-
	if(is.na(le)) " __no length(.)__ "
	else if(give.length) {
	    if(le > 0) paste0("[1:", paste(le), "]") else "(0)"
	} else ""
    v.len <- vec.len # modify v.len, not vec.len!
    ## NON interesting attributes:
    std.attr <- "names"

    cl <- if((S4 <- isS4(object))) class(object) else oldClass(object)
    has.class <- S4 || !is.null(cl) # S3 or S4
    mod <- ""; char.like <- FALSE
    if(give.attr) a <- attributes(object)#-- save for later...
    deParse <- function(.) deparse(., width.cutoff = min(500,max(20, width-10)))

    if (is.null(object))
	cat(" NULL\n")
    else if(S4) {
	if(isRef <- is(object,"envRefClass")) {
	    cld <- object$getClass()
	    nFlds <- names(cld@fieldClasses)
	    a <- sapply(nFlds, function(ch) object[[ch]],
			simplify = FALSE)
	    meths <- ls(cld@refMethods, all.names = TRUE)
	    dfltMs <- ls(getClassDef("envRefClass")@refMethods, all.names = TRUE)
	    oMeths <- meths[is.na(match(meths, dfltMs))]
	    sNms <- names(cld@slots)
	    if(length(sNms <- sNms[sNms != ".xData"]))
		sls <- sapply(sNms, methods::slot,
			      object=object, simplify = FALSE)
	    cat("Reference class", " '", paste(cl, collapse = "', '"),
		"' [package \"", attr(cl,"package"), "\"] with ",
		length(a)," fields\n", sep = "")
	} else {
	    a <- sapply(methods::.slotNames(object), methods::slot,
			object=object, simplify = FALSE)
	    cat("Formal class", " '", paste(cl, collapse = "', '"),
		"' [package \"", attr(cl,"package"), "\"] with ",
		length(a)," slots\n", sep = "")
	}
	if(isRef) {
	    strSub(a, no.list=TRUE, give.length=give.length,
		   nest.lev = nest.lev + 1)
	    cat(indent.str, "and ", length(meths), " methods,", sep = "")
	    if(length(oMeths)) {
		cat(" of which", length(oMeths), "are possibly relevant")
		if (is.na(max.level) || nest.lev < max.level)
		    cat(":",
			strwrap(paste(oMeths, collapse=", "),
				indent = 2, exdent = 2,
				prefix = indent.str, width=width),# exdent = nind),
			sep = "\n")
		else cat("\n")
	    }
	    if(length(sNms)) {
		cat(" and", length(sNms), "slots\n")
		strSub(sls, comp.str = "@ ", no.list=TRUE, give.length=give.length,
		       indent.str = paste(indent.str,".."), nest.lev = nest.lev + 1)
	    }
	}
	else { ## S4
	    strSub(a, comp.str = "@ ", no.list=TRUE, give.length=give.length,
		   indent.str = paste(indent.str,".."), nest.lev = nest.lev + 1)
	}
	return(invisible())
    }
    else if(is.function(object)) {
	cat(if(is.null(ao <- args(object))) deParse(object)
	else { dp <- deParse(ao); paste(dp[-length(dp)], collapse="\n") },"\n")
    } else if(is.list(object)) {
	i.pl <- is.pairlist(object)
	is.d.f <- is.data.frame(object)
	##?if(is.d.f) std.attr <- c(std.attr, "class", if(is.d.f) "row.names")
	if(le == 0) {
	    if(is.d.f) std.attr <- c(std.attr, "class", "row.names")
	    else cat(" ", if(!is.null(names(object))) "Named ",
		     if(i.pl)"pair", "list()\n", sep = "")
	} else { # list, length >= 1 :
	    if(irregCl <- has.class && identical(object[[1L]], object)) {
		le <- length(object <- unclass(object))
		std.attr <- c(std.attr, "class")
	    }
	    if(no.list || (has.class &&
			   any(sapply(paste("str", cl, sep = "."),
					#use sys.function(.) ..
				      function(ob)exists(ob, mode= "function",
							 inherits= TRUE))))) {
		## str.default is a 'NextMethod' : omit the 'List of ..'
		std.attr <- c(std.attr, "class", if(is.d.f) "row.names")
	    } else { # need as.character here for double lengths.
		cat(if(i.pl) "Dotted pair list" else
		    if(irregCl) paste(pClass(cl), "hidden list") else "List",
		    " of ", as.character(le), "\n", sep = "")
	    }
	    if (is.na(max.level) || nest.lev < max.level) {
		nam.ob <-
		    if(is.null(nam.ob <- names(object))) rep.int("", le)
		    else { ncn <- nchar(nam.ob, type="w")
			   if(any(is.na(ncn))) ## slower, but correct:
			      ncn <- vapply(nam.ob, format.info, 0L)
			   format(nam.ob, width = max(ncn), justify="left")
		       }
		for (i in seq_len(min(list.len,le) ) ) {
		    cat(indent.str, comp.str, nam.ob[i], ":", sep = "")
		    envir <- # pass envir for 'promise' components:
			if(typeof(object[[i]]) == "promise") {
			    structure(object, nam= as.name(nam.ob[i]))
			} # else NULL
		    strSub(object[[i]], give.length=give.length,
                           nest.lev = nest.lev + 1,
                           indent.str = paste(indent.str,".."))
		}
	    }
	    if(list.len < le)
		cat(indent.str, "[list output truncated]\n")
	}
    } else { #- not function, not list
	if(is.vector(object)
	   || (is.array(object) && is.atomic(object))
           ## FIXME: is.vector is not documented to allow those modes.
           ## Should this not be is.language?
	   || is.vector(object, mode= "language")
	   || is.vector(object, mode= "symbol")## R bug(<=0.50-a4) should be part
	   ) { ##-- Splus: FALSE for 'named vectors'
	    if(is.atomic(object)) {
		##-- atomic:   numeric	complex	 character  logical
		mod <- substr(mode(object), 1, 4)
		if     (mod == "nume")
		    mod <- if(is.integer(object)) "int"
		    else if(has.class) cl[1L] else "num"
		else if(mod == "char") { mod <- "chr"; char.like <- TRUE }
		else if(mod == "comp") mod <- "cplx" #- else: keep 'logi'
		if(is.array(object)) {
		    rnk <- length(di. <- dim(object))
		    di <- paste0(ifelse(di. > 1, "1:",""), di.,
				 ifelse(di. > 0, "" ," "))
		    pDi <- function(...) paste(c("[", ..., "]"), collapse = "")
		    le.str <- (if(rnk == 1) pDi(di[1L], "(1d)") else
			       pDi(paste0(di[-rnk], ", "), di[rnk]))
		    std.attr <- "dim" #- "names"
		} else if(!is.null(names(object))) {
		    mod <- paste("Named", mod)
		    std.attr <- std.attr[std.attr != "names"]
		}
		if(has.class && length(cl) == 1) {
		    if(cl != mod && substr(cl, 1,nchar(mod)) != mod)
			mod <- paste0("'",cl,"' ", mod)
		    ## don't show the class *twice*
		    std.attr <- c(std.attr, "class")
		}
		str1 <-
		    if(le == 1 && !is.array(object)) paste(NULL, mod)
		    else paste0(" ", mod, if(le>0)" ", le.str)
	    } else { ##-- not atomic, but vector: #
		mod <- typeof(object)#-- typeof(.) is more precise than mode!
		str1 <- switch(mod,
			       call = " call",
			       language = " language",
			       symbol = " symbol",
			       expression = " ",# "expression(..)" by deParse(.)
			       name = " name",
			       ##not in R:argument = "",# .Argument(.) by deParse(.)
			       ## in R (once):	comment.expression

			       ## default :
			       paste("		#>#>", mod, NULL)
			       )
	    }
#  These are S-PLUS classes not found in R.
#	} else if (inherits(object,"rts") || inherits(object,"cts")
#		   || inherits(object,"its")) {
#	    tsp.a <- tspar(object)
#	    t.cl <- cl[b.ts <- substring(cl,2,3) == "ts"] # "rts" "cts" or "its"
#	    ts.kind <- switch(t.cl,
#			      rts="Regular", cts="Calendar", its="Irregular")
#	    ## from  print.summary.ts(.) :
#	    pars <- unlist(sapply(summary(object)$ pars, format,
#				  nsmall=0, digits=digits.d, justify = "none"))
#	    if(length(pars)>=4) pars <- pars[-3]
#	    pars <- paste(abbreviate(names(pars),min=2), pars,
#			  sep= "=", collapse=", ")
#	    str1 <- paste0(ts.kind, " Time-Series ", le.str, " ", pars, ":")
#	    v.len <- switch(t.cl,rts=.8, cts=.6, its=.9) * v.len
#	    class(object) <- if(any(!b.ts)) cl[!b.ts]
#	    std.attr <- c(std.attr, "tspar")
	} else if(stats::is.ts(object)) {
	    tsp.a <- stats::tsp(object)
	    str1 <- paste0(" Time-Series ", le.str, " from ", format(tsp.a[1L]),
			   " to ", format(tsp.a[2L]), ":")
	    std.attr <- c("tsp","class") #- "names"
	} else if (is.factor(object)) {
	    nl <- length(lev.att <- levels(object))
	    if(!is.character(lev.att)) {# should not happen..
		warning("'object' does not have valid levels()")
		nl <- 0
	    } else { ## protect against large nl:
                w <- min(max(width/2, 10), 1000)
                if(nl > w) lev.att <- lev.att[seq_len(w)]
                n.l <- length(lev.att) # possibly  n.l << nl
                lev.att <- encodeString(lev.att, na.encode = FALSE, quote = '"')
            }
	    ord <- is.ordered(object)
	    object <- unclass(object)
	    if(nl) {
		## as from 2.1.0, quotes are included ==> '-2':
		lenl <- cumsum(3 + (nchar(lev.att, type="w") - 2))# level space
		ml <- if(n.l <= 1 || lenl[n.l] <= 13)
		    n.l else which.max(lenl > 13)
		lev.att <- maybe_truncate(lev.att[seq_len(ml)])
	    }
	    else # nl == 0
		ml <- length(lev.att <- "")

	    lsep <- if(ord) "<" else ","
	    str1 <-
		paste0(if(ord)" Ord.f" else " F",
		       "actor w/ ", nl, " level", if(nl != 1) "s",
		       if(nl) " ",
		       if(nl) paste0(lev.att, collapse = lsep),
		       if(ml < nl) paste0(lsep, ".."), ":")

	    std.attr <- c("levels", "class")
	} else if(typeof(object) %in%
		  c("externalptr", "weakref", "environment")) {
	    ## Careful here, we don't want to change pointer objects
	    if(has.class)
                cat(pClass(cl))
	    le <- v.len <- 0
	    str1 <-
		if(is.environment(object)) format(object)
		else paste0("<", typeof(object), ">")
	    has.class <- TRUE # fake for later
	    std.attr <- "class"
	    ## ideally we would figure out if as.character has a
	    ## suitable method and use that.
	} else if(has.class) {
	    cat("Class", if(length(cl) > 1) "es",
		" '", paste(cl, collapse = "', '"), "' ", sep = "")
	    ## If there's a str.<method>, it should have been called before!
	    uo <- unclass(object)
	    if(!is.null(attributes(uo)$class)) {
		## another irregular case
		xtr <- c(if(identical(uo, object)) { # trap infinite loop
		    class(uo) <- NULL
		    "unclass()-immune"
		} else if(!is.object(object)) "not-object")
		if(!is.null(xtr)) cat("{",xtr,"} ", sep = "")
	    }
	    strSub(uo, indent.str = paste(indent.str,".."), nest.lev = nest.lev + 1)
	    return(invisible())
	} else if(is.atomic(object)) {
	    if((1 == length(a <- attributes(object))) && (names(a) == "names"))
		str1 <- paste(" Named vector", le.str)
	    else {
		##-- atomic / not-vector  "unclassified object" ---
		str1 <- paste(" atomic", le.str)
	    }
	} else if(typeof(object) == "promise") {
	    cat(" promise ")
	    if (!is.null(envir)) {
		objExp <- eval(bquote(substitute(.(attr(envir, "nam")), envir)))
		cat("to ")
		strSub(objExp)
	    } else cat(" <...>\n")
	    return(invisible())
	} else {
	    ##-- NOT-atomic / not-vector  "unclassified object" ---
	    ##str1 <- paste(" ??? of length", le, ":")
	    str1 <- paste("length", le)
	}
	##-- end  if else..if else...  {still non-list case}

	##-- This needs some improvement: Not list nor atomic --
	if ((is.language(object) || !is.atomic(object)) && !has.class) {
	    ##-- has.class superfluous --
	    mod <- mode(object)
	    give.mode <- FALSE
	    if (any(mod == c("call", "language", "(", "symbol"))
		|| is.environment(object)) {
		##give.mode <- !is.vector(object)# then it has not yet been done
		if(mod == "(") give.mode <- TRUE
		typ <- typeof(object)
		object <- deParse(object)

		le <- length(object) # is > 1 e.g. for {A;B} language
		format.fun <- function(x)x
		v.len <- round(.5 * v.len)
		if(le > 1 && typ=="language" && object[1L] == "{" && object[le]=="}") {
		    v.len <- v.len + 2
		    if(le >= 3) {
			object <- c(object[1L],
				    paste(sub("^ +", " ", object[2:(le-1)]),
					  collapse = ";"),
				    object[le])
			le <- length(object)
		    }
		}
	    } else if (mod == "expression") {
		format.fun <- function(x) deParse(as.expression(x))
		v.len <- round(.75 * v.len)
	    } else if (mod == "name"){
		object <- paste(object)#-- show `as' char
	    } else if (mod == "argument"){
		format.fun <- deParse
	    } else {
		give.mode <- TRUE
	    }
	    if(give.mode) str1 <- paste0(str1, ', mode "', mod,'":')

	} else if(is.logical(object)) {
	    v.len <- 1.5 * v.len # was '3' originally (but S prints 'T' 'F' ..)
	    format.fun <- formatNum
	} else if(is.numeric(object)) {
	    iv.len <- round(2.5 * v.len)
	    if(iSurv <- inherits(object, "Surv"))
		std.attr <- c(std.attr, "class")
	    int.surv <- iSurv || is.integer(object)
	    if(!int.surv) {
		ob <- if(le > iv.len) object[seq_len(iv.len)] else object
		ao <- abs(ob <- ob[!is.na(ob)])
	    }
	    else if(iSurv)
		le <- length(object <- as.character(object))
	    if(int.surv || (all(ao > 1e-10 | ao==0) && all(ao < 1e10| ao==0) &&
			    all(abs(ob - signif(ob, digits.d)) <= 9e-16*ao))) {
		if(!iSurv || di.[2L] == 2) # "Surv" : implemented as matrix
		    ## use integer-like length
		    v.len <- iv.len
		format.fun <- formatNum
	    } else {
		v.len <- round(1.25 * v.len)
		format.fun <- formatNum
	    }
	} else if(is.complex(object)) {
	    v.len <- round(.75 * v.len)
	    format.fun <- formatNum
	}

	if(char.like) {
	    ## if object is very long, drop the rest which won't be used anyway:
	    max.len <- max(100, width %/% 3 + 1, if(!missing(vec.len)) vec.len)
	    if(le > max.len) object <- object[seq_len(max.len)]
	    encObj <- encodeString(object, quote= '"', na.encode= FALSE)
					#O: encodeString(object)
	    v.len <-
		if(missing(vec.len)) {
		    max(1,sum(cumsum(3 + if(le>0) nchar(encObj, type="w") else 0) <
			      width - (4 + 5*nest.lev + nchar(str1, type="w"))))
		}		      # '5*ne..' above is fudge factor
		else round(v.len)
	    ile <- min(le, v.len)
	    if(ile >= 1) ## truncate if LONG char:
		object <- maybe_truncate(encObj[seq_len(ile)])
					#O: encodeString(object, quote= '"', na.encode= FALSE)
	    formObj <- function(x) paste(as.character(x), collapse=" ")
	}
	else {
	    if(!exists("format.fun", inherits=TRUE)) #-- define one --
		format.fun <-
		    if(mod == "num" || mod == "cplx") format else as.character
	    ## v.len <- max(1,round(v.len))
	    ile <- min(v.len, le)
	    formObj <- function(x) paste(format.fun(x), collapse = " ")
	}

	cat(if(give.head) paste0(str1, " "),
	    formObj(if(ile >= 1) object[seq_len(ile)] else if(v.len > 0) object),
	    if(le > v.len) " ...", "\n", sep = "")

    } ## else (not function nor list)----------------------------------------

    if(give.attr) { ## possible: || has.class && any(cl == "terms")
	nam <- names(a)
	for (i in seq_along(a))
	    if (all(nam[i] != std.attr)) {# only `non-standard' attributes:
		cat(indent.str, paste0('- attr(*, "', nam[i], '")='), sep = "")
		strSub(a[[i]], give.length = give.length,
		       indent.str = paste(indent.str, ".."), nest.lev = nest.lev+1)
	    }
    }
    invisible()	 ## invisible(object)#-- is SLOOOOW on large objects
}# end of `str.default()'

## An extended `ls()' whose print method will use str() :
ls.str <-
    function (pos = -1, name, envir, all.names = FALSE, pattern, mode = "any")
{
    if(missing(envir)) ## [for "lazy" reasons, this fails as default]
        envir <- as.environment(pos)
    nms <- ls(name, envir = envir, all.names=all.names, pattern=pattern)
    r <- unlist(lapply(nms, function(n)
                       exists(n, envir= envir, mode= mode, inherits=FALSE)))
    structure(nms[r], envir = envir, mode = mode, class = "ls_str")
}

lsf.str <- function(pos = -1, envir, ...) {
    if(missing(envir)) ## [for "lazy" reasons, this fails as default]
        envir <- as.environment(pos)
    ls.str(pos=pos, envir=envir, mode = "function", ...)
}

print.ls_str <- function(x, max.level = 1, give.attr = FALSE,
                         ..., digits = max(1, getOption("str")$digits.d))
{
    E <- attr(x, "envir")
    stopifnot(is.environment(E))
    M <- attr(x, "mode")
    args <- list(...)
    if(length(args) && "digits.d" %in% names(args)) {
        if(missing(digits))
            digits <- args$digits.d
        else
            warning("'digits' and 'digits.d' are both specified and the latter is not used")
        args$digits.d <- NULL
    }
    strargs <- c(list(max.level = max.level, give.attr = give.attr,
                      digits = digits), args)
    for(nam in x) {
	cat(nam, ": ")
	## check missingness, e.g. inside debug(.) :

##__ Why does this give	 too many <missing> in some case?
##__	if(eval(substitute(missing(.), list(. = as.name(nam))),
##__		envir = E))
##__	    cat("<missing>\n")
##__	else
##__	    str(get(nam, envir = E, mode = M),
##__		max.level = max.level, give.attr = give.attr, ...)

	o <- tryCatch(get(nam, envir = E, mode = M), error = function(e)e)
	if(inherits(o, "error")) {
	    cat(## FIXME: only works with "C" (or English) LC_MESSAGES locale!
		if(length(grep("missing|not found", o$message)))
		"<missing>" else o$message, "\n", sep = "")
	}
	else {
	    ## do.call(str, c(list(o), strargs),
	    ##	  quote = is.call(o) || is.symbol(o)) # protect calls from eval.
	    ## ensuring 'obj' is *not* copied:
	    strO <- function(...) str(o, ...)
	    do.call(strO, strargs, quote = is.call(o) || is.symbol(o))
					# protect calls from eval.
	}
    }
    invisible(x)
}
#  File src/library/utils/R/summRprof.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

# The profile file always starts with a single header line followed by stack lines
#   If the header contains "memory profiling", the stack lines have memory info
#     The memory info is a fixed width prefix on each line of the form :[0-9]+:[0-9]+:[0-9]+:[0-9]+:
#   If the header contains "line profiling", there will be filename lines and stack lines will contain 
#     line number info of the form [0-9]+#[0-9]+
#   The filename lines will start #File [0-9]+: 

summaryRprof <-
    function(filename = "Rprof.out", chunksize = 5000,
             memory = c("none", "both", "tseries", "stats"),
             lines = c("hide", "show", "both"),
             index = 2, diff = TRUE, exclude = NULL, basenames = 1)
{
    con <- file(filename, "rt")
    on.exit(close(con))
    firstline <- readLines(con, n = 1L)
    if(!length(firstline))
        stop(gettextf("no lines found in %s", sQuote(filename)), domain = NA)
    sample.interval <- as.numeric(strsplit(firstline, "=")[[1L]][2L])/1e6
    memory.profiling <- substr(firstline, 1L, 6L) == "memory"
    line.profiling <- grepl("line profiling", firstline)
    if (line.profiling)
    	filenames <- character(0)

    memory <- match.arg(memory)
    if(memory != "none" && !memory.profiling)
        stop("profile does not contain memory information")
    if (memory == "tseries")
        return(Rprof_memory_summary(filename = con, chunksize = chunksize,
                                    label = index, diff = diff, exclude = exclude,
                                    sample.interval = sample.interval))
    else if (memory == "stats")
        return(Rprof_memory_summary(filename = con,  chunksize = chunksize,
                                    aggregate = index, diff = diff, exclude = exclude,
                                    sample.interval = sample.interval))

    lines <- match.arg(lines)
    if (lines != "hide" && !line.profiling)
    	stop("profile does not contain line information")

    fnames <- NULL
    ucounts <- NULL
    fcounts <- NULL
    memcounts <- NULL
    umem <- NULL

    repeat({

       chunk <- readLines(con, n = chunksize)
 
       if (line.profiling) {
       	   filenamelines <- grep("^#File [0-9]+: ", chunk)
       	   if (length(filenamelines)) {
       	   	fnum <- as.integer(sub("^#File ([0-9]+): .*", "\\1", chunk[filenamelines]))
       	   	filenames[fnum] <- sub("^#File [0-9]+: ", "", chunk[filenamelines])
       	   	if (basenames) {
       		    dirnames <- dirname(filenames[fnum])
       	   	    filenames[fnum] <- basename(filenames[fnum])
       	   	    for (i in seq_len(basenames - 1)) {
       	   	        tail <- basename(dirnames)
       	   	    	filenames[fnum] <- ifelse(tail == ".", filenames[fnum],
       	   	    	                          paste0(tail, "/", filenames[fnum]))
       	   	    	dirnames <- dirname(dirnames)
       	   	    }
       	   	}
       	   	chunk <- chunk[-filenamelines]
       	   }
       }
       
       if (length(chunk) == 0L)
           break
       	       
       if (memory.profiling) {
           memprefix <- attr(regexpr(":[0-9]+:[0-9]+:[0-9]+:[0-9]+:", chunk), "match.length")
           if (memory == "both") {
               memstuff <- substr(chunk, 2L, memprefix-1L)
               memcounts <- pmax(apply(sapply(strsplit(memstuff, ":"), as.numeric), 1, diff), 0)
               ##  memcounts <- c(0, rowSums(memcounts[, 1L:3L]))
               ## convert to bytes.
               memcounts <- c(0, rowSums(cbind(memcounts[, 1L:2L] * 8, memcounts[, 3L])))
               rm(memstuff)
           }
           chunk <- substr(chunk, memprefix+1L, nchar(chunk,  "c"))
           if(any((nc <- nchar(chunk, "c")) == 0L)) {
                chunk <- chunk[nc > 0L]
                memcounts <- memcounts[nc > 0L]
           }
       }

       chunk <- strsplit(chunk, " ")
       if (line.profiling)  
           chunk <- lapply(chunk, function(x) {
           	locations <- !grepl("^\"", x)
           	if (lines != "hide") {
           	    fnum <- sub("#.*", "", x[locations])
           	    lnum <- sub(".*#", "", x[locations])
           	    x[locations] <- paste0(filenames[as.integer(fnum)], "#", lnum)
                }
           	switch(lines,
           	    hide = x <- x[!locations],
           	    show = x <- x[locations]
           	)
       	      	if (length(x)) x else "<no location>"
       	     })
       newfirsts <- sapply(chunk,  "[[",  1L)
       newuniques <- lapply(chunk,  unique)
       ulen <- sapply(newuniques, length)
       newuniques <- unlist(newuniques)

       new.utable <- table(newuniques)
       new.ftable <- table(factor(newfirsts, levels = names(new.utable)))
       if (memory == "both")
           new.umem <- rowsum(memcounts[rep.int(seq_along(memcounts), ulen)], newuniques)

       fcounts <- rowsum( c(as.vector(new.ftable), fcounts),
                         c(names(new.ftable), fnames) )
       ucounts <- rowsum( c(as.vector(new.utable), ucounts),
                         c(names(new.utable), fnames) )
       if(memory == "both")
           umem <- rowsum(c(new.umem, umem), c(names(new.utable), fnames))

       fnames <- sort(unique(c(fnames, names(new.utable))))
    })

    firstnum <- fcounts*sample.interval
    uniquenum <- ucounts*sample.interval

    ## sort and form % on unrounded numbers
    index1 <- order(-firstnum, -uniquenum)
    index2 <- order(-uniquenum, -firstnum)
    
    if (lines == "show") {
    	filename <- sub("#.*$", "", fnames)
    	linenum <- rep(0, length(filename))
    	hasline <- filename != fnames
    	linenum[hasline] <- as.numeric(sub("^.*#", "", fnames[hasline]))
    	index3 <- order(filename, linenum)
    }

    firstpct <- round(100*firstnum/sum(firstnum), 2)
    uniquepct <- round(100*uniquenum/sum(firstnum), 2)

    digits <- ifelse(sample.interval < 0.01,  3L, 2L)
    firstnum <- round(firstnum, digits)
    uniquenum <- round(uniquenum, digits)

    if (memory == "both") memtotal <-  round(umem/1048576, 1)     ## 0.1MB

    rval <- data.frame(firstnum, firstpct, uniquenum, uniquepct)
    names(rval) <- c("self.time", "self.pct", "total.time", "total.pct")
    rownames(rval) <- fnames
    if (memory == "both") rval$mem.total <- memtotal

    by.self <- rval[index1, ]
    by.self <- by.self[by.self[,1L] > 0, ]
    by.total <- rval[index2, c(3L, 4L,  if(memory == "both") 5L, 1L, 2L)]
    
    result <- list(by.self = by.self, by.total = by.total)
    
    if (lines == "show")
    	result <- c(result, list(by.line = rval[index3,]))
    	
    c(result, 
         sample.interval = sample.interval,
         sampling.time = sum(fcounts)*sample.interval)
}

Rprof_memory_summary <- function(filename, chunksize = 5000,
                                 label = c(1, -1), aggregate = 0, diff = FALSE,
                                 exclude = NULL, sample.interval)
{

    fnames <- NULL
    memcounts <- NULL
    firsts <- NULL
    labels <- vector("list", length(label))
    index <- NULL

    repeat({
       chunk <- readLines(filename, n = chunksize)
       if (length(chunk) == 0L)
           break
       memprefix <- attr(regexpr(":[0-9]+:[0-9]+:[0-9]+:[0-9]+:", chunk),
                         "match.length")
       memstuff <- substr(chunk, 2L, memprefix-1L)
       memcounts <- rbind(t(sapply(strsplit(memstuff, ":"), as.numeric)))

       chunk <- substr(chunk, memprefix+1, nchar(chunk,  "c"))
       if(any((nc <- nchar(chunk,  "c")) == 0L)) {
           memcounts <- memcounts[nc > 0L, ]
           chunk <- chunk[nc > 0L]
       }

       chunk <- strsplit(chunk, " ")

       if (length(exclude))
           chunk <- lapply(chunk,  function(l) l[!(l %in% exclude)])

       newfirsts <- sapply(chunk,  "[[",  1L)
       firsts <- c(firsts, newfirsts)

       if (!aggregate && length(label)){
           for(i in seq_along(label)){

               if (label[i] == 1)
                   labels[[i]] <- c(labels[[i]], newfirsts)
               else if (label[i]>1) {
                   labels[[i]] <- c(labels[[i]],  sapply(chunk,
                                                         function(line)
                                                         paste(rev(line)[1L:min(label[i], length(line))],
                                                               collapse = ":")))
               } else {
                   labels[[i]] <- c(labels[[i]], sapply(chunk,
                                                        function(line)
                                                        paste(line[1L:min(-label[i], length(line))],
                                                              collapse = ":")))
               }
           }
       } else if (aggregate) {
           if (aggregate > 0) {
               index <- c(index,  sapply(chunk,
                                         function(line)
                                         paste(rev(line)[1L:min(aggregate, length(line))],
                                               collapse = ":")))

           } else {
               index <- c(index,  sapply(chunk,
                                         function(line)
                                         paste(line[1L:min(-aggregate, length(line))],
                                               collapse = ":")))
           }
       }


       if (length(chunk) < chunksize)
           break
    })

    if (length(memcounts) == 0L) stop("no events were recorded")

    memcounts <- as.data.frame(memcounts)
    names(memcounts) <- c("vsize.small", "vsize.large", "nodes", "duplications")
    if (!aggregate) {
        rownames(memcounts) <- (1L:nrow(memcounts))*sample.interval
        names(labels) <- paste("stack", label, sep = ":")
        memcounts <- cbind(memcounts, labels)
    }

    if (diff)
        memcounts[-1L, 1L:3L]  <-  pmax(0L, apply(memcounts[, 1L:3L], 2L, diff))

    if (aggregate)
        memcounts <- by(memcounts, index,
                        function(these) with(these,
                                             round(c(vsize.small = mean(vsize.small),
                                                     max.vsize.small = max(vsize.small),
                                                     vsize.large = mean(vsize.large),
                                                     max.vsize.large = max(vsize.large),
                                                     nodes = mean(nodes),
                                                     max.nodes = max(nodes),
                                                     duplications = mean(duplications),
                                                     tot.duplications = sum(duplications),
                                                     samples = nrow(these)
                                                     ))
                                             )
                        )
    return(memcounts)
}
#  File src/library/utils/R/tar.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

untar <- function(tarfile, files = NULL, list = FALSE, exdir = ".",
                  compressed = NA, extras = NULL, verbose = FALSE,
                  restore_times = TRUE, tar = Sys.getenv("TAR"))
{
    if (inherits(tarfile, "connection") || identical(tar, "internal"))
        return(untar2(tarfile, files, list, exdir, restore_times))

    if (!(is.character(tarfile) && length(tarfile) == 1L))
        stop("invalid 'tarfile' argument")

    TAR <- tar
    if (!nzchar(TAR) && .Platform$OS.type == "windows" &&
        nzchar(Sys.which("tar.exe"))) TAR <- "tar.exe"
    if (!nzchar(TAR) || TAR == "internal")
        return(untar2(tarfile, files, list, exdir))

    cflag <- ""
    if (is.character(compressed)) {
        ## Any tar which supports -J does not need it for extraction
        switch(match.arg(compressed, c("gzip", "bzip2", "xz")),
               "gzip" = "z", "bzip2" = "j", "xz" = "J")
    } else if (is.logical(compressed)) {
        if (is.na(compressed)) {
            magic <- readBin(tarfile, "raw", n = 3L)
            if(all(magic[1:2] == c(0x1f, 0x8b))) cflag <- "z"
            else if(all(magic[1:2] == c(0x1f, 0x9d))) cflag <- "z" # compress
            else if(rawToChar(magic[1:3]) == "BZh") cflag <- "j"
            else if(rawToChar(magic[1:5]) == "\xFD7zXZ") cflag <- "J"
        } else if (compressed) cflag <- "z"
    } else stop("'compressed' must be logical or character")
    if (!restore_times) cflag <- paste0(cflag, "m")

    gzOK <- .Platform$OS.type == "windows"
    if (!gzOK ) {
        ## version info may be sent to stdout or stderr
        tf <- tempfile()
        ## TAR might be a command+flags, so don't quote it
        cmd <- paste0(TAR, " -", cflag, "tf ", shQuote(tarfile))
        system(paste(TAR, "--version >", tf, "2>&1"))
        if (file.exists(tf)) {
            gzOK <- any(grepl("GNU", readLines(tf), fixed = TRUE))
            unlink(tf)
        }
    }
    tarfile <- path.expand(tarfile)
    if (!gzOK && cflag == "z" && nzchar(ZIP <- Sys.getenv("R_GZIPCMD"))) {
        TAR <- paste(ZIP, "-dc", shQuote(tarfile), "|", TAR)
        tarfile <- "-"
        cflag <- ""
    }
    if (!gzOK && cflag == "j" && nzchar(ZIP <- Sys.getenv("R_BZIPCMD"))) {
        TAR <- paste(ZIP,  "-dc", shQuote(tarfile), "|", TAR)
        tarfile < "-"
        cflag <- ""
    }
    if (cflag == "J") {
        TAR <- paste("xz -dc", shQuote(tarfile), "|", TAR)
        tarfile < "-"
        cflag <- ""
    }
    if (list) {
        cmd <- paste0(TAR, " -", cflag, "tf ", shQuote(tarfile))
        if (length(extras)) cmd <- paste(cmd, extras, collapse = " ")
        if (verbose) message("untar: using cmd = ", sQuote(cmd), domain = NA)
        system(cmd, intern = TRUE)
    } else {
        cmd <- paste0(TAR, " -", cflag, "xf ", shQuote(tarfile))
        if (!missing(exdir)) {
            if (!file_test("-d", exdir)) {
                if(!dir.create(exdir, showWarnings = TRUE, recursive = TRUE))
                    stop(gettextf("failed to create directory %s", sQuote(exdir)),
                         domain = NA)
            }
            cmd <- if(.Platform$OS.type == "windows")
                ## some versions of tar.exe need / here
                paste(cmd, "-C", shQuote(gsub("\\", "/", exdir, fixed=TRUE)))
            else
                paste(cmd, "-C", shQuote(exdir))
        }
        if (length(extras)) cmd <- paste(cmd, extras, collapse = " ")
        if (length(files))
            cmd <- paste(cmd, paste(shQuote(files), collapse = " "))
        if (verbose) message("untar: using cmd = ", sQuote(cmd), domain = NA)
        res <- system(cmd)
        if (res) warning(sQuote(cmd), " returned error code ", res,
                         domain = NA)
        invisible(res)
    }
}

untar2 <- function(tarfile, files = NULL, list = FALSE, exdir = ".",
                   restore_times = TRUE)
{
    ## might be used with len = 12, so result of more than max int
    getOctD <- function(x, offset, len)
    {
        x <- 0.0
        for(i in offset + seq_len(len)) {
            z <- block[i]
            if(!as.integer(z)) break; # terminate on nul
            switch(rawToChar(z),
                   " " = {},
                   "0"=,"1"=,"2"=,"3"=,"4"=,"5"=,"6"=,"7"=
                   {x <- 8*x + (as.integer(z)-48L)},
                   stop("invalid octal digit")
                   )
        }
        x
    }
    getOct <- function(x, offset, len)
        as.integer(getOctD(x, offset, len))
    mydir.create <- function(path, ...) {
        ## for Windows' sake
        path <- sub("[\\/]$", "", path)
        if(file_test("-d", path)) return()
        if(!dir.create(path, showWarnings = TRUE, recursive = TRUE, ...))
           stop(gettextf("failed to create directory %s", sQuote(path)),
                domain = NA)
    }

    warn1 <- character()

    ## A tar file is a set of 512 byte records,
    ## a header record followed by file contents (zero-padded).
    ## See http://en.wikipedia.org/wiki/Tar_%28file_format%29
    if(is.character(tarfile) && length(tarfile) == 1L) {
        con <- gzfile(path.expand(tarfile), "rb") # reads compressed formats
        on.exit(close(con))
    } else if(inherits(tarfile, "connection")) con <- tarfile
    else stop("'tarfile' must be a character string or a connection")
    if (!missing(exdir)) {
        mydir.create(exdir)
        od <- setwd(exdir)
        on.exit(setwd(od), add = TRUE)
    }
    contents <- character()
    llink <- lname <- lsize <- NULL
    repeat{
        block <- readBin(con, "raw", n = 512L)
        if(!length(block)) break
        if(length(block) < 512L) stop("incomplete block on file")
        if(all(block == 0)) break
        ## This should be non-empty, but whole name could be in prefix
        w <- which(block[1:100] > 0)
        ns <- if(length(w)) max(w) else 0
        name <- rawToChar(block[seq_len(ns)])
        magic <- rawToChar(block[258:262])
        if ((magic == "ustar") && block[346L] > 0) {
            ns <- max(which(block[346:500] > 0))
            prefix <- rawToChar(block[345L+seq_len(ns)])
            name <- file.path(prefix, name)
            ns <- nchar(name, "b")
        }
        if (ns <= 0) stop("invalid name field in tarball")
        ## mode zero-padded 8 bytes (including nul) at 101
        ## Aargh: bsdtar has this one incorrectly with 6 bytes+space
        mode <- as.octmode(getOct(block, 100, 8))
        size <- getOctD(block, 124, 12)
        ts <- getOctD(block, 136, 12)
        ft <- as.POSIXct(as.numeric(ts), origin = "1970-01-01", tz = "UTC")
        csum <- getOct(block, 148, 8)
        block[149:156] <- charToRaw(" ")
        xx <- as.integer(block)
        checksum <- sum(xx) %% 2^24 # 6 bytes
        if(csum != checksum) {
            ## try it with signed bytes.
            checksum <- sum(ifelse(xx > 127L, xx - 128L, xx)) %% 2^24 # 6 bytes
            if(csum != checksum)
                warning(gettextf("checksum error for entry '%s'", name),
                        domain = NA)
        }
        type <- block[157L]
        ctype <- rawToChar(type)
#        message(sprintf("%s, %d: '%s'", ctype, size, name))
        if(type %in% c(0L, 7L) || ctype == "0") {
            ## regular or high-performance file
            if(!is.null(lname)) {name <- lname; lname <- NULL}
            if(!is.null(lsize)) {size <- lsize; lsize <- NULL}
            contents <- c(contents, name)
            remain <- size
            dothis <- !list
            if(dothis && length(files)) dothis <- name %in% files
            if(dothis) {
                mydir.create(dirname(name))
                out <- file(name, "wb")
            }
            for(i in seq_len(ceiling(size/512L))) {
                block <- readBin(con, "raw", n = 512L)
                if(length(block) < 512L)
                    stop("incomplete block on file")
                if (dothis) {
                    writeBin(block[seq_len(min(512L, remain))], out)
                    remain <- remain - 512L
                }
            }
            if(dothis) {
                close(out)
                Sys.chmod(name, mode, FALSE) # override umask
                if(restore_times) Sys.setFileTime(name, ft)
            }
        } else if(ctype %in% c("1", "2")) {
            ## hard and symbolic links
            contents <- c(contents, name)
            ns <- max(which(block[158:257] > 0))
            name2 <- rawToChar(block[157L + seq_len(ns)])
            if(!is.null(lname)) {name <- lname; lname <- NULL}
            if(!is.null(llink)) {name2 <- llink; llink <- NULL}
            if(!list) {
                if(ctype == "1") {
                    mydir.create(dirname(name))
                    unlink(name)
                    if (!file.link(name2, name)) { # will give a warning
                        ## link failed, so try a file copy
                        if(file.copy(name2, name))
                             warn1 <- c(warn1, "restoring hard link as a file copy")
                        else
                            warning(gettextf("failed to copy %s to %s", sQuote(name2), sQuote(name)), domain = NA)
                    }
                } else {
                    if(.Platform$OS.type == "windows") {
                        ## this will not work for links to dirs
                        mydir.create(dirname(name))
                        from <- file.path(dirname(name), name2)
                        if (!file.copy(from, name))
                            warning(gettextf("failed to copy %s to %s", sQuote(from), sQuote(name)), domain = NA)
                        else
                            warn1 <- c(warn1, "restoring symbolic link as a file copy")
                   } else {
                       mydir.create(dirname(name))
                       od <- setwd(dirname(name))
                       nm <- basename(name)
                       unlink(nm)
                       if(!file.symlink(name2, nm)) { # will give a warning
                        ## so try a file copy: will not work for links to dirs
                        if (file.copy(name2, nm))
                            warn1 <- c(warn1, "restoring symbolic link as a file copy")
                           else
                               warning(gettextf("failed to copy %s to %s", sQuote(from), sQuote(name)), domain = NA)
                       }
                       setwd(od)
                   }
                }
            }
        } else if(ctype %in% c("3", "4")) {
            ## 3 and 4 are devices
            warn1 <- c(warn1, "skipping devices")
        } else if(ctype == "5") {
            ## directory
            contents <- c(contents, name)
            if(!list) {
                mydir.create(name)
                Sys.chmod(name, mode, TRUE) # FIXME: check result
                ## no point is setting time, as dir will be populated later.
            }
        } else if(ctype == "6") {
            ## 6 is a fifo
            warn1 <- c(warn1, "skipping fifos")
       } else if(ctype %in% c("L", "K")) {
            ## These are GNU extensions that are widely supported
            ## They use one or more blocks to store the name of
            ## a file or link or of a link target.
            name_size <- 512L * ceiling(size/512L)
            block <- readBin(con, "raw", n = name_size)
            if(length(block) < name_size)
                stop("incomplete block on file")
            ns <- max(which(block > 0)) # size on file may or may not include final nul
            if(ctype == "L")
                lname <- rawToChar(block[seq_len(ns)])
            else
                llink <- rawToChar(block[seq_len(ns)])
        } else if(ctype == "x") {
            ## pax headers misused by bsdtar.
            isUTF8 <- FALSE
            warn1 <- c(warn1, "using pax extended headers")
            info <- readBin(con, "raw", n = 512L*ceiling(size/512L))
            info <- strsplit(rawToChar(info), "\n", fixed = TRUE)[[1]]
            hcs <- grep("[0-9]* hdrcharset=", info, useBytes = TRUE,
                        value = TRUE)
            if(length(hcs)) {
                hcs <- sub("[0-9]* hdrcharset=", hcs, useBytes = TRUE)
                isUTF8 <- identical(hcs, "ISO-IR 10646 2000 UTF-8")
            }
            path <- grep("[0-9]* path=", info, useBytes = TRUE, value = TRUE)
            if(length(path)) {
                lname <- sub("[0-9]* path=", "", path, useBytes = TRUE)
                if(isUTF8) Encoding(lname) <- "UTF-8"
            }
            linkpath <- grep("[0-9]* linkpath=", info, useBytes = TRUE,
                             value = TRUE)
            if(length(linkpath)) {
                llink <- sub("[0-9]* linkpath=", "", linkpath, useBytes = TRUE)
                if(isUTF8) Encoding(llink) <- "UTF-8"
            }
            size <- grep("[0-9]* size=", info, useBytes = TRUE, value = TRUE)
            if(length(size))
                lsize <- as.integer(sub("[0-9]* size=", "", size))
         } else if(ctype == "g") {
            warn1 <- c(warn1, "skipping pax global extended headers")
            readBin(con, "raw", n = 512L*ceiling(size/512L))
        } else stop("unsupported entry type ", sQuote(ctype))
    }
    if(length(warn1)) {
        warn1 <- unique(warn1)
        for (w in warn1) warning(w, domain = NA)
    }
    if(list) contents else invisible(0L)
}

tar <- function(tarfile, files = NULL,
                compression = c("none", "gzip", "bzip2", "xz"),
                compression_level = 6, tar = Sys.getenv("tar"),
                extra_flags = "")
{
    if(is.character(tarfile)) {
        if(nzchar(tar) && tar != "internal") {
            ## FIXME: could pipe through gzip etc: might be safer for xz
            ## as -J was lzma in GNU tar 1.20:21
            flags <- switch(match.arg(compression),
                            "none" = "-cf",
                            "gzip" = "-zcf",
                            "bzip2" = "-jcf",
                            "xz" = "-Jcf")

            if (grepl("darwin", R.version$os)) {
                ## precaution for Mac OS X to omit resource forks
                ## we can't tell the running OS version from R.version$os
                ## but at least it will not be older
                tar <- paste("COPYFILE_DISABLE=1", tar) # >= 10.5, Leopard
                if (grepl("darwin8", R.version$os)) # 10.4, Tiger
                    tar <- paste("COPY_EXTENDED_ATTRIBUTES_DISABLE=1", tar)
            }
            if (is.null(extra_flags)) extra_flags <- ""
            ## 'tar' might be a command + flags, so don't quote it
            cmd <- paste(tar, extra_flags, flags, shQuote(tarfile),
                         paste(shQuote(files), collapse=" "))
            return(invisible(system(cmd)))
        }
        con <- switch(match.arg(compression),
                      "none" =    file(tarfile, "wb"),
                      "gzip" =  gzfile(tarfile, "wb", compression = compression_level),
                      "bzip2" = bzfile(tarfile, "wb", compression = compression_level),
                      "xz" =    xzfile(tarfile, "wb", compression = compression_level))
        on.exit(close(con))
    } else if(inherits(tarfile, "connection")) con <- tarfile
    else stop("'tarfile' must be a character string or a connection")

    ## FIXME: eventually we should use the pax extension, but
    ## that was first supported in R 2.15.3.
    GNUname <- function(name, link = FALSE)
    {
        header <- raw(512L)
        n1 <- charToRaw("ExtendedName")
        header[seq_along(n1)] <- n1
        header[157L] <- charToRaw(ifelse(link, "K", "L"))
        size <- length(name)
        header[125:135] <- charToRaw(sprintf("%011o", as.integer(size)))
        header[149:156] <- charToRaw(" ")
        checksum <- sum(as.integer(header)) %% 2^24 # 6 bytes
        header[149:154] <- charToRaw(sprintf("%06o", as.integer(checksum)))
        header[155L] <- as.raw(0L)
        writeBin(header, con)
        writeBin(name, con)
        ssize <- 512L * ceiling(size/512L)
        if(ssize > size) writeBin(raw(ssize - size), con)
    }
    warn1 <- character()

    files <- list.files(files, recursive = TRUE, all.files = TRUE,
                        full.names = TRUE, include.dirs = TRUE)

    for (f in unique(files)) {
        info <- file.info(f)
        if(is.na(info$size)) {
            warning(gettextf("file '%s' not found", f), domain = NA)
            next
        }
        header <- raw(512L)
        ## add trailing / to dirs.
        if(info$isdir && !grepl("/$", f)) f <- paste0(f, "/")
        name <- charToRaw(f)
        if(length(name) > 100L) {
            OK <- TRUE
            ## best possible case: 155+/+100
            if(length(name) > 256L) OK <- FALSE
            else {
                ## do not want to split on terminal /
                m <- length(name)
                s <- max(which(name[1:min(156, m - 1L)] == charToRaw("/")))
                if(is.infinite(s) || s + 100L < length(name)) OK <- FALSE
            }
            warning("storing paths of more than 100 bytes is not portable:\n  ",
                    sQuote(f), domain = NA)
            if (OK) {
                prefix <- name[1:(s-1L)]
                name <- name[-(1:s)]
                header[345L+seq_along(prefix)] <- prefix
            } else {
                GNUname(name)
                name <- charToRaw("dummy")
                warn1 <- c(warn1, "using GNU extension for long pathname")
            }
        }
        header[seq_along(name)] <- name
        mode <- info$mode
        ## for use by R CMD build
        if (is.null(extra_flags) && grepl("/(configure|cleanup)$", f) &&
            (mode & "111") != as.octmode("111")) {
            warning(gettextf("file '%s' did not have execute permissions: corrected", f), domain = NA, call. = FALSE)
            mode <- mode | "111"
        }
        header[101:107] <- charToRaw(sprintf("%07o", mode))
        ## Windows does not have uid, gid: defaults to 0, which isn't great
        uid <- info$uid
        ## uids are supposed to be less than 'nobody' (32767)
        ## but it seems there are broken ones around: PR#15436
        if(!is.null(uid) && !is.na(uid)) {
            if(uid < 0L || uid > 32767L) {
                warning(gettextf("invalid uid value replaced by that for user 'nobody'", uid),
                        domain = NA, call. = FALSE)
                uid <- 32767L
            }
            header[109:115] <- charToRaw(sprintf("%07o", uid))
        }
        gid <- info$gid
        if(!is.null(gid) && !is.na(gid)) {
            if(gid < 0L || gid > 32767L) {
                warning(gettextf("invalid gid value replaced by that for user 'nobody'", uid),
                        domain = NA, call. = FALSE)
                gid <- 32767L
            }
            header[117:123] <- charToRaw(sprintf("%07o", gid))
	}
        header[137:147] <- charToRaw(sprintf("%011o", as.integer(info$mtime)))
        if (info$isdir) header[157L] <- charToRaw("5")
        else {
            lnk <- Sys.readlink(f)
            if(is.na(lnk)) lnk <- ""
            header[157L] <- charToRaw(ifelse(nzchar(lnk), "2", "0"))
            if(nzchar(lnk)) {
                if(nchar(lnk, "b") > 100L) {
                    ##  stop("linked path is too long")
                    GNUname(charToRaw(lnk), TRUE)
                    warn1 <- c(warn1, "using GNU extension for long linkname")
                    lnk <- "dummy"
                }
                header[157L + seq_len(nchar(lnk))] <- charToRaw(lnk)
                size <- 0
            }
        }
        ## size is 0 for directories and it seems for links.
        size <- ifelse(info$isdir, 0, info$size)
        if(size >= 8^11) stop("file size is limited to 8GB")
        header[125:135] <- .Call(C_octsize, size)
        ## the next two are what POSIX says, not what GNU tar does.
        header[258:262] <- charToRaw("ustar")
        header[264:265] <- charToRaw("0")
        ## Windows does not have uname, grname
        s <- info$uname
        if(!is.null(s) && !is.na(s)) {
            ns <- nchar(s, "b")
            header[265L + (1:ns)] <- charToRaw(s)
        }
        s <- info$grname
        if(!is.null(s) && !is.na(s)) {
            ns <- nchar(s, "b")
            header[297L + (1:ns)] <- charToRaw(s)
        }
        header[149:156] <- charToRaw(" ")
        checksum <- sum(as.integer(header)) %% 2^24 # 6 bytes
        header[149:154] <- charToRaw(sprintf("%06o", as.integer(checksum)))
        header[155L] <- as.raw(0L)
        writeBin(header, con)
        if(info$isdir || nzchar(lnk)) next
        inf <- file(f, "rb")
        for(i in seq_len(ceiling(info$size/512L))) {
            block <- readBin(inf, "raw", 512L)
            writeBin(block, con)
            if( (n <- length(block)) < 512L) writeBin(raw(512L - n), con)
        }
        close(inf)
    }
    ## trailer is two blocks of nuls.
    block <- raw(512L)
    writeBin(block, con)
    writeBin(block, con)
    if(length(warn1)) {
        warn1 <- unique(warn1)
        for (w in warn1) warning(w, domain = NA)
    }
    invisible(0L)
}
#  File src/library/utils/R/toLatex.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

toBibtex <- function(object, ...) UseMethod("toBibtex")

toLatex <- function(object, ...) UseMethod("toLatex")

print.Bibtex <- print.Latex <- function(x, prefix = "", ...)
{
    writeLines(paste0(prefix, unclass(x)), ...)
    invisible(x)
}
#  File src/library/utils/R/unix/create.post.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

create.post <- function(instructions = character(),
                        description = "post",
                        subject = "",
                        method = getOption("mailer"),
                        address = "the relevant mailing list",
                        ccaddress = getOption("ccaddress", ""),
                        filename = "R.post",
                        info = character())
{
    method <-
	if(is.null(method)) "none"
	else match.arg(method, c("mailto", "mailx", "gnudoit", "none", "ess"))
    open_prog <- if(grepl("-apple-darwin", R.version$platform)) "open" else "xdg-open"
    if (method == "mailto")
        if(!nzchar(Sys.which(open_prog))) {
            browser <- Sys.getenv("R_BROWSER", "")
            if(!nzchar(browser)) {
                warning("cannot find program to open 'mailto:' URIs: reverting to 'method=\"none\"'")
                flush.console()
                Sys.sleep(5)
            } else {
                message("Using the browser to open a mailto: URI")
                open_prog <- browser
            }
        }

    body <- c(instructions,
              "--please do not edit the information below--", "",
              info)

    none_method <- function() {
        disclaimer <-
            paste0("# Your mailer is set to \"none\",\n",
                   "# hence we cannot send the, ", description, " directly from R.\n",
                   "# Please copy the ", description, " (after finishing it) to\n",
                   "# your favorite email program and send it to\n#\n",
                   "#       ", address, "\n#\n",
                   "######################################################\n",
                   "\n\n")

        cat(c(disclaimer, body), file = filename, sep = "\n")
        cat("The", description, "is being opened for you to edit.\n")
        flush.console()
        file.edit(filename)
        cat("The unsent ", description, " can be found in file ",
            sQuote(filename), "\n", sep ="")
    }

    if(method == "none") {
        none_method()
    } else if(method == "mailx") {
        if(missing(address)) stop("must specify 'address'")
        if(!nzchar(subject)) stop("'subject' is missing")
        if(length(ccaddress) != 1L) stop("'ccaddress' must be of length 1")

	cat(body, file=filename, sep = "\n")
        cat("The", description, "is being opened for you to edit.\n")
        file.edit(filename)

        if(is.character(ccaddress) && nzchar(ccaddress)) {
            cmdargs <- paste("-s", shQuote(subject),
                             "-c", shQuote(ccaddress),
                             shQuote(address),
                             "<", filename, "2>/dev/null")
        }
        else
            cmdargs <- paste("-s", shQuote(subject),
                             shQuote(address), "<",
                             filename, "2>/dev/null")
        status <- 1L
        answer <- readline(paste0("Email the ", description, " now? (yes/no) "))
        answer <- grep("yes", answer, ignore.case=TRUE)
        if(length(answer)) {
            cat("Sending email ...\n")
            status <- system(paste("mailx", cmdargs), , TRUE, TRUE)
            if(status)
                status <- system(paste("Mail", cmdargs), , TRUE, TRUE)
            if(status)
                status <- system(paste("/usr/ucb/mail", cmdargs), , TRUE, TRUE)

            if(status == 0L) unlink(filename)
            else {
                cat("Sending email failed!\n")
                cat("The unsent", description, "can be found in file",
                    sQuote(filename), "\n")
            }
        } else
            cat("The unsent", description, "can be found in file", filename, "\n")
    } else if(method == "ess") {
	cat(body, sep = "\n")
    } else  if(method == "gnudoit") {
        ## FIXME: insert subject and ccaddress
	cmd <- paste0("gnudoit -q '",
		     "(mail nil \"", address, "\")",
		     "(insert \"", paste(body, collapse="\\n"), "\")",
		     "(search-backward \"Subject:\")",
		     "(end-of-line)'")
	system(cmd)
    } else if(method == "mailto") {
        if (missing(address)) stop("must specify 'address'")
        if (!nzchar(subject)) subject <- "<<Enter Meaningful Subject>>"
        if(length(ccaddress) != 1L) stop("'ccaddress' must be of length 1")
        cat("The", description, "is being opened in your default mail program\nfor you to complete and send.\n")
        ## The mailto: standard (RFC2368) says \r\n for the body
        arg <- paste0("mailto:", address,
                     "?subject=", subject,
                     if(is.character(ccaddress) && nzchar(ccaddress))
                         paste0("&cc=", ccaddress),
                     "&body=", paste(body, collapse = "\r\n"))
        if(system2(open_prog, shQuote(URLencode(arg)), FALSE, FALSE)) {
            cat("opening the mailer failed, so reverting to 'mailer=\"none\"'\n")
            flush.console()
            none_method()
        }
    }
    invisible()
}
#  File src/library/utils/R/unix/download.file.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

download.file <-
    function(url, destfile, method, quiet = FALSE, mode = "w",
             cacheOK = TRUE, extra = getOption("download.file.extra"))
{
    destfile # check supplied
    method <- if (missing(method))
	getOption("download.file.method", default = "auto")
    else
        match.arg(method, c("auto", "internal", "wget", "curl", "lynx"))

    if(method == "auto") {
        if(capabilities("http/ftp"))
            method <- "internal"
        else if(length(grep("^file:", url))) {
            method <- "internal"
            url <- URLdecode(url)
        } else if(system("wget --help > /dev/null") == 0L)
            method <- "wget"
        else if(system("curl --help > /dev/null") == 0L)
            method <- "curl"
        else if(system("lynx -help > /dev/null") == 0L)
            method <- "lynx"
        else
            stop("no download method found")
    }
    if(method == "internal") {
        status <- .External(C_download, url, destfile, quiet, mode, cacheOK)
        ## needed for Mac GUI from download.packages etc
        if(!quiet) flush.console()
    } else if(method == "wget") {
        if(quiet) extra <- c(extra, "--quiet")
        if(!cacheOK) extra <- c(extra, "--cache=off")
        status <- system(paste("wget",
                               paste(extra, collapse = " "),
                               shQuote(url),
                               "-O", shQuote(path.expand(destfile))))
    } else if(method == "curl") {
        if(quiet) extra <- c(extra, "-s -S")
        if(!cacheOK) extra <- c(extra, "-H 'Pragma: no-cache'")
        status <- system(paste("curl",
                               paste(extra, collapse = " "),
                               shQuote(url),
                               " -o", shQuote(path.expand(destfile))))
    } else if(method == "lynx")
        status <- system(paste("lynx -dump",
                               paste(extra, collapse = " "),
                               shQuote(url),
                               ">", shQuote(path.expand(destfile))))

    if(status) warning("download had nonzero exit status")

    invisible(status)
}

nsl <- function(hostname) .Call(C_nsl, hostname)
#  File src/library/utils/R/unix/mac.install.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


if(substr(R.version$os, 1L, 6L) != "darwin") {
.install.macbinary <-
    function(pkgs, lib, repos = getOption("repos"),
             contriburl = contrib.url(repos, type="mac.binary"),
             method, available = NULL, destdir = NULL,
             dependencies = FALSE,
             lock = getOption("install.lock", FALSE), quiet = FALSE,
             ...)
    {}
} else {
## edited from windows/.install.winbinary
##
.install.macbinary <-
    function(pkgs, lib, repos = getOption("repos"),
             contriburl = contrib.url(repos, type="mac.binary"),
             method, available = NULL, destdir = NULL,
             dependencies = FALSE,
             lock = getOption("install.lock", FALSE), quiet = FALSE,
             ...)
{
    untar <- function(what, where)
    {
        ## FIXME: should this look for Sys.getenv('TAR')?
        ## Leopard has GNU tar, SL has BSD tar.
        xcode <- system(paste0("tar zxf \"", path.expand(what), "\" -C \"",
                               path.expand(where), "\""), intern=FALSE)
        if (xcode)
            warning(gettextf("'tar' returned non-zero exit code %d", xcode),
                    domain = NA, call. = FALSE)
    }

    unpackPkg <- function(pkg, pkgname, lib, lock = FALSE)
    {
        dir.exists <- function(x) !is.na(isdir <- file.info(x)$isdir) & isdir
        ## Create a temporary directory and unpack the zip to it
        ## then get the real package & version name, copying the
        ## dir over to the appropriate install dir.
        tmpDir <- tempfile(, lib)
        if (!dir.create(tmpDir))
            stop(gettextf("unable to create temporary directory %s",
                          sQuote(tmpDir)),
                 domain = NA, call. = FALSE)
        cDir <- getwd()
        on.exit(setwd(cDir), add = TRUE)
        res <- untar(pkg, tmpDir)
        setwd(tmpDir)
        ## sanity check: people have tried to install source .tgz files
        if (!file.exists(file <- file.path(pkgname, "Meta", "package.rds")))
            stop(gettextf("file %s is not an OS X binary package", sQuote(pkg)),
                 domain = NA, call. = FALSE)
        desc <- readRDS(file)$DESCRIPTION
        if (length(desc) < 1L)
            stop(gettextf("file %s is not an OS X binary package", sQuote(pkg)),
                 domain = NA, call. = FALSE)
        desc <- as.list(desc)
        if (is.null(desc$Built))
            stop(gettextf("file %s is not an OS X binary package", sQuote(pkg)),
                 domain = NA, call. = FALSE)

        res <- tools::checkMD5sums(pkgname, file.path(tmpDir, pkgname))
        if(!quiet && !is.na(res) && res) {
            cat(gettextf("package %s successfully unpacked and MD5 sums checked\n",
                         sQuote(pkgname)))
            flush.console()
        }

        instPath <- file.path(lib, pkgname)
        if(identical(lock, "pkglock") || isTRUE(lock)) {
	    lockdir <- if(identical(lock, "pkglock"))
                file.path(lib, paste("00LOCK", pkgname, sep = "-"))
            else file.path(lib, "00LOCK")
	    if (file.exists(lockdir)) {
                stop(gettextf("ERROR: failed to lock directory %s for modifying\nTry removing %s",
                              sQuote(lib), sQuote(lockdir)), domain = NA)
	    }
	    dir.create(lockdir, recursive = TRUE)
	    if (!dir.exists(lockdir))
                stop(gettextf("ERROR: failed to create lock directory %s",
                              sQuote(lockdir)), domain = NA)
            ## Back up a previous version
            if (file.exists(instPath)) {
                file.copy(instPath, lockdir, recursive = TRUE)
        	on.exit({
         	    if (restorePrevious) {
                        try(unlink(instPath, recursive = TRUE))
        	    	savedcopy <- file.path(lockdir, pkgname)
        	    	file.copy(savedcopy, lib, recursive = TRUE)
        	    	warning(gettextf("restored %s", sQuote(pkgname)),
                                domain = NA, call. = FALSE, immediate. = TRUE)
        	    }
        	}, add=TRUE)
        	restorePrevious <- FALSE
            }
	    on.exit(unlink(lockdir, recursive = TRUE), add=TRUE)
        }
        ## If the package is already installed, remove it.  If it
        ## isn't there, the unlink call will still return success.
        ret <- unlink(instPath, recursive=TRUE)
        if (ret == 0L) {
            ## Move the new package to the install lib and
            ## remove our temp dir
            ret <- file.rename(file.path(tmpDir, pkgname), instPath)
            if(!ret) {
                warning(gettextf("unable to move temporary installation %s to %s",
                                 sQuote(file.path(tmpDir, pkgname)),
                                 sQuote(instPath)),
                        domain = NA, call. = FALSE)
                restorePrevious <- TRUE # Might not be used
            }
        } else
        stop(gettextf("cannot remove prior installation of package %s",
                      sQuote(pkgname)), call. = FALSE, domain = NA)
        setwd(cDir)
        unlink(tmpDir, recursive=TRUE)
    }

    if(!length(pkgs)) return(invisible())

    if(is.null(contriburl)) {
        pkgnames <- basename(pkgs)
        pkgnames <- sub("\\.tgz$", "", pkgnames)
        pkgnames <- sub("\\.tar\\.gz$", "", pkgnames)
        pkgnames <- sub("_.*$", "", pkgnames)
        ## there is no guarantee we have got the package name right:
        ## foo.zip might contain package bar or Foo or FOO or ....
        ## but we can't tell without trying to unpack it.
        for(i in seq_along(pkgs))
            unpackPkg(pkgs[i], pkgnames[i], lib, lock = lock)
        return(invisible())
    }
    tmpd <- destdir
    nonlocalcran <- length(grep("^file:", contriburl)) < length(contriburl)
    if(is.null(destdir) && nonlocalcran) {
        tmpd <- file.path(tempdir(), "downloaded_packages")
        if (!file.exists(tmpd) && !dir.create(tmpd))
            stop(gettextf("unable to create temporary directory %s",
                          sQuote(tmpd)),
                 domain = NA)
    }

    if(is.null(available))
        available <- available.packages(contriburl = contriburl,
                                        method = method)
    pkgs <- getDependencies(pkgs, dependencies, available, lib)

    foundpkgs <- download.packages(pkgs, destdir = tmpd, available = available,
                                   contriburl = contriburl, method = method,
                                   type = "mac.binary", quiet = quiet, ...)

    if(length(foundpkgs)) {
        update <- unique(cbind(pkgs, lib))
        colnames(update) <- c("Package", "LibPath")
        for(lib in unique(update[,"LibPath"])) {
            oklib <- lib==update[,"LibPath"]
            for(p in update[oklib, "Package"])
            {
                okp <- p == foundpkgs[, 1L]
                if(any(okp))
                    unpackPkg(foundpkgs[okp, 2L], foundpkgs[okp, 1L], lib,
                              lock = lock)
            }
        }
        if(!quiet && !is.null(tmpd) && is.null(destdir))
            cat("\n", gettextf("The downloaded binary packages are in\n\t%s", tmpd),
                "\n", sep = "")
    } else if(!is.null(tmpd) && is.null(destdir)) unlink(tmpd, recursive = TRUE)

    invisible()
}
}
#  File src/library/utils/R/unix/sysutils.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

memory.size <- function(max = FALSE)
{
    warning("'memory.size()' is Windows-specific", call.=FALSE)
    Inf
}

memory.limit <- function(size = NA)
{
   warning("'memory.limit()' is Windows-specific", call.=FALSE)
   Inf
}
#  File src/library/utils/R/utils-defunct.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

## <entry>
## Deprecated in 1.9.0
## Defunct in 2.0.0
## Removed in 3.0.0
## package.contents <- function(pkg, lib.loc=NULL) .Defunct(package="utils")
## </entry>

## <entry>
## Deprecated in 2.12.2
## Defunct in 2.14.0
## Removed in 3.0.0
## zip.file.extract <- function(file, zipname = "R.zip",
## 			     unzip = getOption("unzip"), dir = tempdir())
## .Defunct("unzip")
## </entry>

## <entry>
## Deprecated in 2.2.0
## Defunct in 3.0.0
CRAN.packages <- function(CRAN = getOption("repos"), method,
                          contriburl = contrib.url(CRAN))
    .Defunct("available.packages")
## </entry>
#  File src/library/utils/R/utils-deprecated.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#  File src/library/utils/R/vignette.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

vignette <-
    function(topic, package = NULL, lib.loc = NULL, all = TRUE)
{
    vinfo <- tools:::getVignetteInfo(package, lib.loc, all)
    
    if(!missing(topic)) {
        topic <- topic[1L]               # Just making sure ...
        vinfo <- vinfo[vinfo[, "Topic"] == topic,,drop=FALSE]
        if(length(vinfo)) {

            pdf <- vinfo[, "PDF"]
            pidx <- file_test("-f", file.path(vinfo[, "Dir"], "doc", vinfo[, "PDF"]))

            if(any(pidx)){
                idx <- min(which(pidx))
                if(sum(pidx)>1){
                    ## <FIXME>
                    ## Should really offer a menu to select from.
                    warning(gettextf("vignette %s found more than once,\nusing the one found in %s",
                                     sQuote(topic), sQuote(dirname(pdf[idx]))),
                            call. = FALSE, domain = NA)
                    ## </FIXME>
                }
		vinfo <- vinfo[idx,,drop=FALSE]
		Dir <- vinfo[, "Dir"]
		File <- vinfo[, "File"]
		PDF <- vinfo[, "PDF"]
                z <- list(file=file.path(Dir, "doc", File),
                          pdf=file.path(Dir, "doc", PDF))
            }
            else{
		Dir <- vinfo[1, "Dir"]
		File <- vinfo[1, "File"]
                z <- list(file=file.path(Dir, "doc", File),
                          pdf=character(0L))
            }
            z$topic <- topic
            class(z) <- "vignette"
            return(z)
        }
        else
            warning(gettextf("vignette %s not found", sQuote(topic)),
                    call. = FALSE, domain = NA)
    }

    if(missing(topic)) {
        ## List all possible vignettes.

        title <- if(nrow(vinfo)) {
            paste(vinfo[, "Title"],
                  paste0(rep.int("(source", nrow(vinfo)),
                        ifelse(vinfo[, "PDF"] != "", paste0(", ", tools::file_ext(vinfo[, "PDF"])), ""),
                        ")"))
        }
        else
            character()
        ## ... and rewrite into the form used by packageIQR.
        db <- cbind(Package = basename(vinfo[, "Dir"]),
                    LibPath = dirname(vinfo[, "Dir"]),
                    Item = vinfo[, "Topic"],
                    Title = title)
	footer <- if (all) NULL else
		  paste0("Use ",
                         sQuote("vignette(all = TRUE)"),
                         "\n",
                         "to list the vignettes in all *available* packages.")

        y <- list(type = "vignette", title = "Vignettes", header = NULL,
                  results = db, footer = footer)
        class(y) <- "packageIQR"
        return(y)
    }
}

print.vignette <- function(x, ...){

    if(length(x$pdf)){
        ## <FIXME>
        ## Should really abstract this into a BioC style
        ## openPDF() along the lines of browseURL() ...
        ext <- tools::file_ext(x$pdf)
        if (tolower(ext) == "pdf") {
            pdfviewer <- getOption("pdfviewer")
            if(identical(pdfviewer, "false")) {
            } else if(.Platform$OS.type == "windows" &&
                      identical(pdfviewer, file.path(R.home("bin"), "open.exe")))
            	shell.exec(x$pdf)
            else system2(pdfviewer, shQuote(x$pdf), wait = FALSE)
        ## </FIXME>         
        } else 
             browseURL(x$pdf)

    } else {
        warning(gettextf("vignette %s has no PDF/HTML", sQuote(x$topic)),
                call. = FALSE, domain = NA)
    }
    invisible(x)
}

edit.vignette <- function(name, ...)
{

    f <- tempfile(name$topic, fileext=".R")
    Stangle(name$file, output=f, quiet=TRUE)
    file.edit(file=f, ...)
}
#  File src/library/utils/R/widgets.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

select.list <-
    function(choices, preselect = NULL, multiple = FALSE, title = NULL,
             graphics = getOption("menu.graphics"))
{
    if(!interactive()) stop("select.list() cannot be used non-interactively")
    if(!is.null(title) && (!is.character(title) || length(title) != 1))
        stop("'title' must be NULL or a length-1 character vector")
    if(isTRUE(graphics)) {
        if (.Platform$OS.type == "windows" || .Platform$GUI == "AQUA")
        return(.External2(C_selectlist, choices, preselect, multiple, title))
        ## must be Unix here
        ## Tk might not require X11 on Mac OS X, but if DISPLAY is set
        ## this will work for Aqua Tcl/Tk.
        ## OTOH, we do want to check Tk works!
        else if(graphics && capabilities("tcltk") &&
                capabilities("X11") && suppressWarnings(tcltk:::.TkUp))
            return(tcltk::tk_select.list(choices, preselect, multiple, title))
    }
    ## simple text-based alternatives.
    if(!multiple) {
        res <- menu(choices, FALSE, title)
        if(res < 1L || res > length(choices)) return("")
        else return(choices[res])
    } else {
        nc <- length(choices)
        if (length(title) && nzchar(title[1L]))
            cat(title, "\n", sep = "")
        def <- if(is.null(preselect)) rep(FALSE, nc)
        else choices %in% preselect
        op <- paste0(format(seq_len(nc)), ": ",
                     ifelse(def, "+", " "), " ", choices)
        if(nc > 10L) {
            fop <- format(op)
            nw <- nchar(fop[1L], "w") + 2L
            ncol <- getOption("width") %/% nw
            if(ncol > 1L)
                op <- paste(fop, c(rep("  ", ncol - 1L), "\n"),
                            sep = "", collapse="")
            cat("", op, sep = "\n")
        } else cat("", op, "", sep = "\n")
        cat(gettext("Enter one or more numbers separated by spaces, or an empty line to cancel\n"))
	repeat {
            res <- tryCatch(scan("", what = 0, quiet = TRUE, nlines = 1),
                            error = identity)
	    if(!inherits(res, "error")) break
	    cat(gettext("Invalid input, please try again\n"))
	}
        if(!length(res) || (length(res) == 1L && !res[1L])) return(character())
        res <- sort(res[1 <= res && res <= nc])
        return(choices[res])
    }
}

flush.console <- function() invisible(.Call(C_flushconsole))

process.events <- function() invisible(.Call(C_processevents))
#  File src/library/utils/R/write.table.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

write.table <-
function (x, file = "", append = FALSE, quote = TRUE, sep = " ",
          eol = "\n", na = "NA", dec = ".", row.names = TRUE,
          col.names = TRUE, qmethod = c("escape", "double"),
          fileEncoding = "")
{
    qmethod <- match.arg(qmethod)
    if(is.logical(quote) && (length(quote) != 1L || is.na(quote)))
        stop("'quote' must be 'TRUE', 'FALSE' or numeric")
    ## quote column names unless quote == FALSE (see help).
    quoteC <- if(is.logical(quote)) quote else TRUE
    qset <- is.logical(quote) && quote

    if(!is.data.frame(x) && !is.matrix(x)) x <- data.frame(x)

    makeRownames <- isTRUE(row.names)
    ## need col names if col.names is TRUE or NA
    makeColnames <- is.logical(col.names) && !identical(FALSE, col.names)
    if(is.matrix(x)) {
        ## fix up dimnames as as.data.frame would
        p <- ncol(x)
        d <- dimnames(x)
        if(is.null(d)) d <- list(NULL, NULL)
        if(is.null(d[[1L]]) && makeRownames) d[[1L]] <- seq_len(nrow(x))
        if(is.null(d[[2L]]) && makeColnames && p > 0L)
            d[[2L]] <- paste0("V", 1L:p)
        if(qset)
            quote <- if(is.character(x)) seq_len(p) else numeric()
    } else { ## data.frame
        if(qset)
            quote <- if(length(x))
                which(unlist(lapply(x, function(x)
                                    is.character(x) || is.factor(x))))
            else numeric()
        ## fix up embedded matrix columns into separate cols:
        if(any(sapply(x, function(z) length(dim(z)) == 2 && dim(z)[2L] > 1))) {
            c1 <- names(x)
	    x <- as.matrix(x, rownames.force = makeRownames)
	    d <- dimnames(x)
	    if(qset) {
		ord <- match(c1, d[[2L]], 0L)
		quote <- ord[quote]; quote <- quote[quote > 0L]
	    }
        }
        else
            d <- list(if(makeRownames) row.names(x),
                      if(makeColnames) names(x))
        p <- ncol(x)
    }
    nocols <- p == 0L

    if(is.logical(quote)) # must be false
	quote <- NULL
    else if(is.numeric(quote)) {
	if(any(quote < 1L | quote > p))
	    stop("invalid numbers in 'quote'")
    } else
	stop("invalid 'quote' specification")

    rn <- FALSE
    rnames <- NULL
    if(is.logical(row.names)) {
	if(row.names) {rnames <- as.character(d[[1L]]); rn <- TRUE}
    } else {
	rnames <- as.character(row.names)
        rn <- TRUE
	if(length(rnames) != nrow(x))
            stop("invalid 'row.names' specification")
    }
    if(!is.null(quote) && rn) # quote the row names
	quote <- c(0, quote)

    if(is.logical(col.names)) {
        if(!rn && is.na(col.names))
            stop("'col.names = NA' makes no sense when 'row.names = FALSE'")
        col.names <- if(is.na(col.names) && rn) c("", d[[2L]])
        else if(col.names) d[[2L]] else NULL
    } else {
	col.names <- as.character(col.names)
	if(length(col.names) != p)
	    stop("invalid 'col.names' specification")
    }

    if(file == "") file <- stdout()
    else if(is.character(file)) {
        file <- if(nzchar(fileEncoding))
            file(file, ifelse(append, "a", "w"), encoding = fileEncoding)
            else file(file, ifelse(append, "a", "w"))
        on.exit(close(file))
    } else if(!isOpen(file, "w")) {
        open(file, "w")
        on.exit(close(file))
    }
    if(!inherits(file, "connection"))
        stop("'file' must be a character string or connection")

    qstring <-                          # quoted embedded quote string
        switch(qmethod,
               "escape" = '\\\\"',
               "double" = '""')
    if(!is.null(col.names)) {
	if(append)
	    warning("appending column names to file")
	if(quoteC)
	    col.names <- paste("\"", gsub('"', qstring, col.names),
                               "\"", sep = "")
        writeLines(paste(col.names, collapse = sep), file, sep = eol)
    }

    if (nrow(x) == 0L) return(invisible())
    if (nocols && !rn) return(cat(rep.int(eol, NROW(x)), file=file, sep=""))

    ## convert list matrices to character - maybe not much use?
    if(is.matrix(x) && !is.atomic(x)) mode(x) <- "character"
    if(is.data.frame(x)) {
        ## convert columns we can't handle in C code
        x[] <- lapply(x, function(z) {
            if(is.object(z) && !is.factor(z)) as.character(z) else z
        })
    }

    invisible(.External2(C_writetable, x, file, nrow(x), p, rnames, sep, eol,
                         na, dec, as.integer(quote), qmethod != "double"))
}

write.csv <- function(...)
{
    Call <- match.call(expand.dots = TRUE)
    for(argname in c("append", "col.names", "sep", "dec", "qmethod"))
        if(!is.null(Call[[argname]]))
            warning(gettextf("attempt to set '%s' ignored", argname),
                    domain = NA)
    rn <- eval.parent(Call$row.names)
    Call$append <- NULL
    Call$col.names <- if(is.logical(rn) && !rn) TRUE else NA
    Call$sep <- ","
    Call$dec <- "."
    Call$qmethod <- "double"
    Call[[1L]] <- as.name("write.table")
    eval.parent(Call)
}

write.csv2 <- function(...)
{
    Call <- match.call(expand.dots = TRUE)
    for(argname in c("append", "col.names", "sep", "dec", "qmethod"))
        if(!is.null(Call[[argname]]))
            warning(gettextf("attempt to set '%s' ignored", argname),
                    domain = NA)
    rn <- eval.parent(Call$row.names)
    Call$append <- NULL
    Call$col.names <- if(is.logical(rn) && !rn) TRUE else NA
    Call$sep <- ";"
    Call$dec <- ","
    Call$qmethod <- "double"
    Call[[1L]] <- as.name("write.table")
    eval.parent(Call)
}
#  File src/library/utils/R/zip.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

unzip <-
    function(zipfile, files = NULL, list = FALSE, overwrite = TRUE,
             junkpaths = FALSE, exdir = ".", unzip = "internal",
             setTimes = FALSE)
{
    if(identical(unzip, "internal")) {
        if(!list && !missing(exdir))
            dir.create(exdir, showWarnings = FALSE, recursive = TRUE)
        res <- .External(C_unzip, zipfile, files, exdir, list, overwrite,
                         junkpaths, setTimes)
        if(list) {
            dates <- as.POSIXct(res[[3]], "%Y-%m-%d %H:%M",  tz="UTC")
            data.frame(Name = res[[1]], Length = res[[2]], Date = dates,
                       stringsAsFactors = FALSE)
        } else invisible(attr(res, "extracted"))
    } else {
        WINDOWS <- .Platform$OS.type == "windows"
        if(!is.character(unzip) || length(unzip) != 1L || !nzchar(unzip))
            stop("'unzip' must be a single character string")
        zipfile <- path.expand(zipfile)
        if (list) {
            res <- if (WINDOWS)
                    system2(unzip, c("-l", shQuote(zipfile)), stdout = TRUE)
                else
                    system2(unzip, c("-l", shQuote(zipfile)), stdout = TRUE,
                            env = c("TZ=UTC"))
            l <- length(res)
            res2 <- res[-c(1,3, l-1, l)]
            con <- textConnection(res2); on.exit(close(con))
            z <- read.table(con, header=TRUE, as.is=TRUE)
            dt <- paste(z$Date, z$Time)
            ## Unzip 6.00 always uses 4-digits years, but any order is
            ## possible and the separator could be - or / (depending
            ## on the locale on Windows).
            ## Unzip 5.52 uses 2-digit years, but default to "%m-%d-%y" on
            ## most platforms (but is locale-dependent on Windows).
            formats <-
                if (max(nchar(z$Date) > 8))
                    c("%Y-%m-%d", "%d-%m-%Y", "%m-%d-%Y") else
                    ## At this point we are guessing: there is no way
                    ## to know what "08-09-10" means.  Take the most common
                    ## default first.
                    c("%m-%d-%y", "%d-%m-%y", "%y-%m-%d")
            slash <- any(grepl("/", z$Date))
            if (slash) formats <- gsub("-", "/", formats)
            formats <- paste(formats, "%H:%M")
            for (f in formats) {
                zz <- as.POSIXct(dt, tz="UTC", format = f)
                if (all(!is.na(zz))) break
            }
            z[, "Date"] <- zz
            z[c("Name", "Length", "Date")]
        } else {
            args <- c("-oq", shQuote(zipfile))
            if (length(files)) args <- c(args, shQuote(files))
            if (exdir != ".") args <- c(args, "-d", shQuote(exdir))
            ## there is an unzip clone about that does not respect -q
            system2(unzip, args, stdout = NULL, stderr = NULL,
                    invisible = TRUE)
            invisible(NULL)
        }
    }
}

zip <- function(zipfile, files, flags = "-r9X", extras = "",
                zip = Sys.getenv("R_ZIPCMD", "zip"))
{
    if (missing(flags) && (!is.character(files) || !length(files)))
        stop("'files' must a character vector specifying one or more filepaths")
    args <- c(flags, shQuote(path.expand(zipfile)),
              shQuote(files), extras)
    if (.Platform$OS.type == "windows")
        invisible(system2(zip, args, invisible = TRUE))
    else invisible(system2(zip, args))
}

#  File src/library/utils/R/zzz.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

.noGenerics <- TRUE

.onLoad <- function(libname, pkgname)
{
    ## Set default options() related to functionality in 'utils' pkg
    op <- options()
    op.utils <-
	list(help.try.all.packages = FALSE,
	     help.search.types = c("vignette", "demo", "help"),
             citation.bibtex.max = 1, internet.info = 2,
	     pkgType = .Platform$pkgType,
	     str = list(strict.width = "no", digits.d = 3, vec.len = 4),
	     demo.ask = "default", example.ask = "default",
	     HTTPUserAgent = defaultUserAgent(),
	     menu.graphics = TRUE, mailer = "mailto")
    extra <-
        if(.Platform$OS.type == "windows") {
            list(unzip = "internal",
                 editor = if(length(grep("Rgui", commandArgs(), TRUE))) "internal" else "notepad",
                 repos = c(CRAN="@CRAN@",
                           CRANextra="http://www.stats.ox.ac.uk/pub/RWin")
                 )
        } else
            list(unzip = Sys.getenv("R_UNZIPCMD"),
                 editor = Sys.getenv("EDITOR"),
                 repos = c(CRAN="@CRAN@"))
    op.utils <- c(op.utils, extra)
    toset <- !(names(op.utils) %in% names(op))
    if(any(toset)) options(op.utils[toset])
}
