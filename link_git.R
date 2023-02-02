install.packages("usethis")
usethis::use_git_config(user.name= "HeleneGSK",
                        user.email="helene.x.fourmanoir@gsk.com")
usethis::use_git()
usethis::create_github_token()



usethis::use_github()
