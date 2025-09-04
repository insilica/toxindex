import os
import uuid
from pathlib import Path

import psycopg2

from webserver.storage import GCSFileStorage
from webserver.model.task import Task
from workflows.celery_task_ranking import ranking_task


ASSISTANT_USER_ID = "00000000-0000-0000-0000-000000000000"


def get_db_connection():
    return psycopg2.connect(
        host=os.getenv("PGHOST"),
        port=os.getenv("PGPORT"),
        dbname=os.getenv("PGDATABASE"),
        user=os.getenv("PGUSER"),
        password=os.getenv("PGPASSWORD"),
    )


def ensure_workflow(title: str) -> int:
    with get_db_connection() as con:
        con.autocommit = True
        with con.cursor() as cur:
            # Try to find existing workflow first
            cur.execute("SELECT workflow_id FROM workflows WHERE title = %s", (title,))
            row = cur.fetchone()
            if row:
                return int(row[0])

            # Fix sequence misalignment if needed (common after manual imports)
            cur.execute(
                """
                SELECT setval(
                  pg_get_serial_sequence('workflows','workflow_id'),
                  COALESCE((SELECT MAX(workflow_id) FROM workflows), 0),
                  true
                )
                """
            )

            # Insert new row and return id
            cur.execute(
                "INSERT INTO workflows (title) VALUES (%s) RETURNING workflow_id",
                (title,),
            )
            return int(cur.fetchone()[0])


def insert_file(file_id: str, task_id: str, user_id: str, filename: str, gcs_path: str) -> None:
    with get_db_connection() as con:
        con.autocommit = True
        with con.cursor() as cur:
            cur.execute(
                """
                INSERT INTO files (file_id, task_id, user_id, filename, filepath)
                VALUES (%s, %s, %s, %s, %s)
                """,
                (file_id, task_id, user_id, filename, gcs_path),
            )


def main():
    workflow_id = ensure_workflow("ranking")
    print("workflow_id:", workflow_id)

    task = Task.create_task(
        title="Ranking Test",
        user_id=ASSISTANT_USER_ID,
        workflow_id=workflow_id,
        description="GHS nephrotoxicity ranking test",
    )
    if not task:
        raise RuntimeError("Task.create_task failed")
    print("task_id:", task.task_id)

    local_path = "/app/ranking-workflow/ranking_workflow/data/pdaa_phthalate_inchis_sample.txt"
    filename = Path(local_path).name
    gcs_path = f"uploads/{task.task_id}/{filename}"
    GCSFileStorage().upload_file(local_path, gcs_path, content_type="text/plain")
    file_id = str(uuid.uuid4())
    insert_file(file_id, str(task.task_id), ASSISTANT_USER_ID, filename, gcs_path)
    print("input_file_id:", file_id)

    payload = {
        "task_id": str(task.task_id),
        "user_id": ASSISTANT_USER_ID,
        "payload": file_id,
        "score_type": "GHS",
        "endpoint": "nephrotoxicity",
    }
    res = ranking_task.delay(payload)
    print("celery_id:", res.id)


if __name__ == "__main__":
    main()


