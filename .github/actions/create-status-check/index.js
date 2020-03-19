const core = require("@actions/core");
const github = require("@actions/github");

async function run() {
  try {
    const [
      gitHubRepoOwner,
      gitHubRepoName
    ] = process.env.GITHUB_REPOSITORY.split("/");
    const gitHubSha = process.env.GITHUB_SHA;
    const gitHubToken = core.getInput("github-token");

    const checkName   = core.getInput("name");
    const checkStatus = core.getInput("status"); # 'in_progress', 'completed'
    const checkState  = core.getInput("state");

    check = {
      owner: gitHubRepoOwner,
      repo: gitHubRepoName,
      name: checkName,
      head_sha: gitHubSha,
      status: checkStatus,
      conclusion: checkState,
    };

    if (checkStatus == 'in_progress') {
        check = {
            ..check,
            ..{
                started_at: new Date().toISOString()
            }
        };
    }

    const outputTitle   = core.getInput("output-title");
    const outputSummary = core.getInput("output-summary");
    if (outputTitle) {
        check = {
            ..check,
            ..{
                title: outputTitle,
                summary: outputSummary
            }
        };
    }

    const octokit = new github.GitHub(gitHubToken);

    octokit.checks.create(check);

    core.setOutput("time", new Date().toTimeString());
  } catch (error) {
    core.setFailed(error.message);
  }
}

run();
