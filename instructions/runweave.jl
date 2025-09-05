using Weave
kwargs=(doctype = "md2html", template = "math2504assessment.tpl")
weave("project1.jmd"; kwargs...) 

