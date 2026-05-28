import os
import sys
import json
import subprocess
import urllib.request
import urllib.error

def run_git_cmd(args):
    """Execute a git command and return stdout/stderr."""
    try:
        res = subprocess.run(args, capture_output=True, text=True, check=True)
        return res.stdout.strip(), ""
    except subprocess.CalledProcessError as e:
        return "", e.stderr.strip()

def make_github_api_call(url, data=None, token=None, method="GET"):
    """Perform a GitHub REST API call using urllib."""
    req = urllib.request.Request(url, method=method)
    req.add_header("Accept", "application/vnd.github+json")
    req.add_header("User-Agent", "RDKit-Studio-Setup")
    if token:
        req.add_header("Authorization", f"Bearer {token}")
    
    if data:
        json_data = json.dumps(data).encode("utf-8")
        req.add_header("Content-Type", "application/json")
        req.data = json_data
        
    try:
        with urllib.request.urlopen(req) as response:
            res_body = response.read().decode("utf-8")
            status = response.status
            return status, json.loads(res_body) if res_body else {}
    except urllib.error.HTTPError as e:
        err_body = e.read().decode("utf-8")
        try:
            err_json = json.loads(err_body)
        except Exception:
            err_json = {"message": err_body}
        return e.code, err_json
    except Exception as e:
        return 500, {"message": str(e)}

def main():
    print("==================================================")
    print("🧪 RDKit Studio: GitHub Push & Collaborator Inviter")
    print("==================================================")
    print("This script will help you push this project to GitHub and add a collaborator.")
    print("Please make sure you have generated a Personal Access Token (PAT) with 'repo' scope.")
    print("Create a PAT at: https://github.com/settings/tokens")
    print("--------------------------------------------------")
    
    # Inputs
    username = input("Enter your GitHub username: ").strip()
    if not username:
        print("Error: GitHub username is required.")
        sys.exit(1)
        
    token = input("Enter your Personal Access Token (PAT): ").strip()
    if not token:
        print("Error: Personal Access Token is required.")
        sys.exit(1)
        
    collab = input("Enter collaborator's GitHub username to invite (e.g. shekar gudda): ").strip()
    if not collab:
        print("Error: Collaborator username is required.")
        sys.exit(1)
        
    repo_name = "rdkit-studio"
    print("\n[1/4] Creating repository on GitHub...")
    
    # 1. Create Repository
    create_url = "https://api.github.com/user/repos"
    create_data = {
        "name": repo_name,
        "description": "Interactive Chemical Intelligence Dashboard using Streamlit and RDKit",
        "private": False
    }
    
    status, res = make_github_api_call(create_url, data=create_data, token=token, method="POST")
    
    if status == 201:
        print(f"✅ Repository '{repo_name}' created successfully on GitHub!")
    elif status == 422 and "name already exists" in str(res.get("errors", "")):
        print(f"ℹ️ Repository '{repo_name}' already exists on your GitHub. Using existing one.")
    else:
        print(f"❌ Failed to create repository (Status {status}): {res.get('message', 'Unknown error')}")
        if status == 401:
            print("Please double check that your Personal Access Token (PAT) is correct and has 'repo' scopes.")
        sys.exit(1)
        
    print("\n[2/4] Setting up local Git remotes...")
    # 2. Add Remote
    # Remove existing remote if any
    run_git_cmd(["git", "remote", "remove", "origin"])
    
    # Add authenticated remote
    remote_url = f"https://{username}:{token}@github.com/{username}/{repo_name}.git"
    stdout, stderr = run_git_cmd(["git", "remote", "add", "origin", remote_url])
    if stderr:
        print(f"❌ Failed to configure remote: {stderr}")
        sys.exit(1)
    print("✅ Configured remote origin successfully.")
    
    print("\n[3/4] Pushing code to GitHub (main branch)...")
    # 3. Push code
    # Streamlit and packages are excluded by .gitignore
    stdout, stderr = run_git_cmd(["git", "push", "-u", "origin", "main"])
    if "Failed" in stderr or "error" in stderr:
        # Sometimes git prints progress to stderr, so we verify by return code in run_git_cmd which raises if push fails
        print(f"❌ Failed to push code: {stderr}")
        sys.exit(1)
    else:
        print("✅ Code successfully pushed to GitHub branch 'main'!")
        
    print("\n[4/4] Sending collaborator invitation...")
    # 4. Invite Collaborator
    collab_url = f"https://api.github.com/repos/{username}/{repo_name}/collaborators/{collab}"
    status, res = make_github_api_call(collab_url, token=token, method="PUT")
    
    if status == 201:
        print(f"🎉 Successfully invited '{collab}' to collaborate on the repository!")
        print("An invitation email has been sent, and they can accept it at their GitHub notifications.")
    elif status == 204:
        print(f"ℹ️ '{collab}' is already a collaborator on this repository.")
    else:
        print(f"❌ Failed to invite collaborator (Status {status}): {res.get('message', 'Unknown error')}")
        
    print("\n==================================================")
    print("🎉 GitHub setup complete!")
    print(f"Your repository is live at: https://github.com/{username}/{repo_name}")
    print("==================================================")

if __name__ == "__main__":
    main()
