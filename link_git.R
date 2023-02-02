install.packages("usethis")
usethis::use_git_config(user.name= "hf69578",
                        user.email="helene.x.fourmanoir@gsk.com")
usethis::use_git()
usethis::create_github_token()



gitcreds::gitcreds_set()

usethis::use_github()
