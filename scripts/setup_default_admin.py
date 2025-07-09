#!/usr/bin/env python3
"""
Script to set up a default admin user.
This should be run after the Flyway migrations have created the schema.
"""

import os
import sys
import webserver.datastore as ds
from werkzeug.security import generate_password_hash

def setup_default_admin(email='test@test.com', password='test'):
    """Set up default admin user."""
    # Check if admin user already exists
    existing_user = ds.find("SELECT user_id FROM users WHERE email = %s", (email,))
    
    if existing_user:
        print(f"Admin user {email} already exists")
        return
    
    # Check if admin group exists
    admin_group = ds.find("SELECT group_id FROM user_groups WHERE name = 'admin'")
    
    if not admin_group:
        print("Admin group not found. Please run migrations first.")
        return
    
    admin_group_id = admin_group['group_id']
    
    # Create admin user
    import uuid
    user_id = uuid.uuid4()
    hashpw = generate_password_hash(password)
    
    ds.execute("""
        INSERT INTO users (user_id, email, hashpw, token, stripe_customer_id, email_verified, group_id)
        VALUES (%s, %s, %s, %s, %s, %s, %s)
    """, (
        user_id,
        email,
        hashpw,
        'default_admin_token',
        'default_admin_stripe',
        True,  # Email verified
        admin_group_id
    ))
    
    print(f"Default admin user created:")
    print(f"  Email: {email}")
    print(f"  Password: {password}")
    print(f"  User ID: {user_id}")
    print(f"  Group: admin")

def main():
    """Main function."""
    try:
        # Set up default admin
        setup_default_admin()
        print("Default admin setup completed successfully")
        
    except Exception as e:
        print(f"Database error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 