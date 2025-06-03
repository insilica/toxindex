import webserver.datastore as ds
import logging


class Environment:
    def __init__(
        self, environment_id, title, user_id, description=None, created_at=None
    ):
        self.environment_id = environment_id
        self.title = title
        self.user_id = user_id
        self.description = description
        self.created_at = created_at

    def to_dict(self):
        return {
            "environment_id": self.environment_id,
            "title": self.title,
            "user_id": self.user_id,
            "description": self.description,
            "created_at": (
                self.created_at.strftime("%Y-%m-%d %H:%M:%S")
                if self.created_at
                else None
            ),
        }

    @staticmethod
    def from_row(row):
        return Environment(
            environment_id=row["environment_id"],
            title=row["title"],
            user_id=row["user_id"],
            description=row.get("description"),
            created_at=row["created_at"],
        )

    @staticmethod
    def create_environment(title, user_id, description=None):
        params = (title, user_id, description)
        ds.execute(
            "INSERT INTO environments (title, user_id, description) VALUES (%s, %s, %s)",
            params,
        )
        logging.info(f"created environment {title} for user {user_id}")
        res = ds.find(
            "SELECT * FROM environments WHERE title = %s AND user_id = %s ORDER BY created_at DESC LIMIT 1",
            (title, user_id),
        )
        return Environment.from_row(res) if res else None

    @staticmethod
    def get_environments_by_user(user_id):
        rows = ds.find_all(
            "SELECT * FROM environments WHERE user_id = %s ORDER BY created_at DESC",
            (user_id,),
        )
        return [Environment.from_row(row) for row in rows]

    @staticmethod
    def get_environment(environment_id):
        row = ds.find(
            "SELECT * FROM environments WHERE environment_id = %s", (environment_id,)
        )
        return Environment.from_row(row) if row else None

    @staticmethod
    def delete_environment(environment_id, user_id=None):
        params = (environment_id,)
        query = "DELETE FROM environments WHERE environment_id = %s"
        if user_id is not None:
            query += " AND user_id = %s"
            params = (environment_id, user_id)
        ds.execute(query, params)
        logging.info(f"Deleted environment {environment_id} for user {user_id}")
