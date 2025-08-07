#!/usr/bin/env python3
"""
Test script to check workflow database state and identify the issue.
"""

import os
import sys
import webserver.datastore as ds
from webserver.model.user_group import UserGroup
from webserver.model.workflow import Workflow

def test_database_connection():
    """Test if database connection is working."""
    print("Testing database connection...")
    try:
        # Test basic connection
        result = ds.find("SELECT 1 as test")
        if result:
            print("✓ Database connection successful")
            return True
        else:
            print("✗ Database connection failed")
            return False
    except Exception as e:
        print(f"✗ Database connection error: {e}")
        return False

def test_workflows_in_database():
    """Check if workflows exist in the database."""
    print("\nTesting workflows in database...")
    try:
        workflows = ds.find_all("SELECT * FROM workflows")
        print(f"Found {len(workflows)} workflows in database:")
        for w in workflows:
            print(f"  - ID: {w['workflow_id']}, Title: {w['title']}")
        return workflows
    except Exception as e:
        print(f"✗ Error querying workflows: {e}")
        return []

def test_user_groups():
    """Check user groups in database."""
    print("\nTesting user groups...")
    try:
        groups = UserGroup.get_all_groups()
        print(f"Found {len(groups)} user groups:")
        for g in groups:
            print(f"  - ID: {g.group_id}, Name: {g.name}")
        return groups
    except Exception as e:
        print(f"✗ Error querying user groups: {e}")
        return []

def test_workflow_access():
    """Check workflow access permissions."""
    print("\nTesting workflow access permissions...")
    try:
        access_records = ds.find_all("SELECT * FROM workflow_group_access")
        print(f"Found {len(access_records)} workflow access records:")
        for a in access_records:
            print(f"  - Group: {a['group_id']}, Workflow: {a['workflow_id']}")
        return access_records
    except Exception as e:
        print(f"✗ Error querying workflow access: {e}")
        return []

def test_get_accessible_workflows(user_id):
    """Test the get_accessible_workflows function."""
    print(f"\nTesting get_accessible_workflows for user {user_id}...")
    try:
        accessible = UserGroup.get_accessible_workflows(user_id)
        print(f"Found {len(accessible)} accessible workflows:")
        for w in accessible:
            print(f"  - ID: {w['workflow_id']}, Title: {w['title']}")
        return accessible
    except Exception as e:
        print(f"✗ Error getting accessible workflows: {e}")
        return []

def main():
    """Main test function."""
    print("=== Workflow Database Test ===")
    
    # Test database connection
    if not test_database_connection():
        print("Cannot proceed without database connection")
        return
    
    # Test workflows
    workflows = test_workflows_in_database()
    
    # Test user groups
    groups = test_user_groups()
    
    # Test workflow access
    access_records = test_workflow_access()
    
    # Test accessible workflows for a sample user (you might need to adjust this)
    if groups:
        # Try with the first user group
        sample_user_id = "00000000-0000-0000-0000-000000000001"  # admin user
        test_get_accessible_workflows(sample_user_id)
    
    print("\n=== Test Summary ===")
    print(f"Database connection: {'✓' if test_database_connection() else '✗'}")
    print(f"Workflows in DB: {len(workflows)}")
    print(f"User groups: {len(groups)}")
    print(f"Access records: {len(access_records)}")
    
    if not workflows:
        print("\n⚠️  No workflows found in database. You may need to run:")
        print("   python scripts/seed_workflows.py")
    
    if not access_records:
        print("\n⚠️  No workflow access records found. You may need to run:")
        print("   python scripts/seed_workflows.py")

if __name__ == "__main__":
    main()
