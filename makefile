git_push_all:
	git push && git switch release && git merge main && git push && git switch main


test_latest:
	./scripts/test_versions.bash latest
