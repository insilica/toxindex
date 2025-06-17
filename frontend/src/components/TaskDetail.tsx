import React, { useEffect, useState, useCallback } from "react";
import { useParams, useNavigate } from "react-router-dom";

interface Task {
  task_id: string;
  title: string;
  archived?: boolean;
  environment_id?: string;
  created_at?: string;
  finished_at?: string;
  session_id?: string;
  user_id?: string;
  workflow_id?: number;
  // Add other fields as needed
}

const getDuration = (start?: string, end?: string) => {
  if (!start || !end) return null;
  const startTime = new Date(start).getTime();
  const endTime = new Date(end).getTime();
  const diff = Math.max(0, endTime - startTime);
  const mins = Math.floor(diff / 60000);
  const secs = Math.floor((diff % 60000) / 1000);
  return `${mins}m ${secs}s`;
};

const TaskDetail: React.FC = () => {
  const { task_id } = useParams<{ task_id: string }>();
  const [task, setTask] = useState<Task | null>(null);
  const [loading, setLoading] = useState(true);
  const navigate = useNavigate();

  const fetchTask = useCallback(async () => {
    if (!task_id) return;
    setLoading(true);
    try {
      const res = await fetch(`/api/tasks/${task_id}`, { credentials: "include" });
      if (!res.ok) {
        setTask(null);
        setLoading(false);
        return;
      }
      const data = await res.json();
      if (data.error) {
        setTask(null);
      } else {
        setTask(data);
      }
    } catch (e) {
      setTask(null);
    } finally {
      setLoading(false);
    }
  }, [task_id]);

  useEffect(() => {
    fetchTask();
  }, [fetchTask]);

  if (loading) return <div className="text-white p-8">Loading...</div>;
  if (!task) return <div className="text-white p-8">Task not found.</div>;

  const duration = getDuration(task.created_at, task.finished_at);

  return (
    <div className="max-w-xl mx-auto p-8 bg-gray-900 rounded-lg shadow text-white mt-12">
      <h2 className="text-2xl font-bold mb-4">{task.title}</h2>
      <div className="mb-2"><b>Tool:</b> {task.workflow_id === 1 ? "ToxIndex RAP" : task.workflow_id === 2 ? "ToxIndex Vanilla" : task.workflow_id ?? <span className="text-gray-400">Unknown</span>}</div>
      <div className="mb-2"><b>Executed by:</b> {task.user_id ? (
        <button
          className="text-purple-400 underline hover:text-purple-300 focus:outline-none"
          onClick={() => navigate(`/user/${task.user_id}`)}
          title="Go to user profile"
          aria-label="Go to user profile"
        >
          {task.user_id}
        </button>
      ) : (
        <span className="text-gray-400">Unknown</span>
      )}</div>
      <div className="mb-2"><b>Created:</b> {task.created_at ? new Date(task.created_at).toLocaleString() : <span className="text-gray-400">Unknown</span>}</div>
      {task.finished_at && <div className="mb-2"><b>Finished:</b> {new Date(task.finished_at).toLocaleString()}</div>}
      {duration && <div className="mb-2"><b>Duration:</b> {duration}</div>}
      <div className="mb-2"><b>Environment:</b> {task.environment_id ? (
        <button
          className="text-blue-400 underline hover:text-blue-300 focus:outline-none"
          onClick={() => navigate(`/environment/${task.environment_id}`)}
          title="Go to environment"
          aria-label="Go to environment"
        >
          {task.environment_id}
        </button>
      ) : (
        <span className="text-gray-400">None</span>
      )}</div>
      <div className="mb-2"><b>Chat Session:</b> {task.session_id ? (
        <button
          className="text-green-400 underline hover:text-green-300 focus:outline-none"
          onClick={() => navigate(`/chat/session/${task.session_id}`)}
          title="Go to chat session"
          aria-label="Go to chat session"
        >
          {task.session_id}
        </button>
      ) : (
        <span className="text-gray-400">None</span>
      )}</div>
      <button className="mt-4 px-4 py-2 bg-green-700 rounded" onClick={() => navigate(-1)} aria-label="Back">Back</button>
    </div>
  );
};

export default TaskDetail; 